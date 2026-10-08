path.distance <- function(path){
dists=c()
for(i in 2:nrow(path)){
x1 = path[i-1,1]
y1 = path[i-1,2]
x2 = path[i,1]
y2 = path[i,2]
d.x = x2 - x1
d.y = y2 - y1
dist = sqrt(d.x*d.x + d.y*d.y)
dists = append(dist,dists)
}
return(mean(dists))
}

densify_lineage_graph <- function(cds, lineage, spacing = NULL, factor = 1.5,
                                  reduction_method = "UMAP", reselect_cells = FALSE,
                                  update_principal_graph = FALSE,
                                  N = 5, cl = 1, sel_clusters = NULL,
                                  start_regions = FALSE, starting_clusters = FALSE) {
  
  if (is.null(cds@graphs[[lineage]]))
    stop("No graph for lineage '", lineage, "' in cds@graphs.", call. = FALSE)
  
  g_sub <- cds@graphs[[lineage]]
  Y     <- cds@principal_graph_aux[[reduction_method]]$dp_mst      # dims x nodes
  vs    <- V(g_sub)$name
  miss  <- setdiff(vs, colnames(Y))
  if (length(miss) > 0)
    stop("Lineage nodes missing from dp_mst: ", paste(miss, collapse = ", "), call. = FALSE)
  
  coords <- t(Y[, vs, drop = FALSE])                              # nodes x dims
  el     <- get.edgelist(g_sub)                                   # E x 2 (names)
  if (nrow(el) == 0) stop("Lineage graph '", lineage, "' has no edges.", call. = FALSE)
  
  edge_len <- sqrt((coords[el[, 1], 1] - coords[el[, 2], 1])^2 +
                     (coords[el[, 1], 2] - coords[el[, 2], 2])^2)
  if (is.null(spacing)) spacing <- stats::median(edge_len)
  thresh <- factor * spacing
  
  # next free Y_<k> index across the WHOLE dp_mst (keeps new names globally unique)
  all_idx <- suppressWarnings(as.integer(sub("^Y_", "", colnames(Y))))
  next_id <- max(all_idx[is.finite(all_idx)], 0L) + 1L
  
  new_names   <- character(0)
  new_coords  <- matrix(numeric(0), nrow = 0, ncol = 2)
  new_edge_df <- data.frame(from = character(0), to = character(0), stringsAsFactors = FALSE)
  
  for (e in seq_len(nrow(el))) {
    a <- el[e, 1]; b <- el[e, 2]; L <- edge_len[e]
    n_seg <- if (L > thresh) max(2L, as.integer(round(L / spacing))) else 1L
    if (n_seg == 1L) {
      new_edge_df <- rbind(new_edge_df, data.frame(from = a, to = b, stringsAsFactors = FALSE))
      next
    }
    A <- coords[a, ]; B <- coords[b, ]
    fracs <- seq_len(n_seg - 1L) / n_seg
    ins   <- paste0("Y_", next_id + seq_len(n_seg - 1L) - 1L)
    next_id <- next_id + (n_seg - 1L)
    xy <- cbind(A[1] + fracs * (B[1] - A[1]),
                A[2] + fracs * (B[2] - A[2]))
    rownames(xy) <- ins
    new_names   <- c(new_names, ins)
    new_coords  <- rbind(new_coords, xy)
    chain       <- c(a, ins, b)
    new_edge_df <- rbind(new_edge_df,
                         data.frame(from = utils::head(chain, -1),
                                    to   = utils::tail(chain, -1),
                                    stringsAsFactors = FALSE))
  }
  
  n_subdiv <- sum(edge_len > thresh)
  if (length(new_names) == 0) {
    message(sprintf("No edges exceeded %.3g (= %.3g x spacing %.3g); nothing to densify.",
                    thresh, factor, spacing))
    return(cds)
  }
  message(sprintf("Densified lineage '%s': added %d node(s) across %d edge(s) (spacing %.3g).",
                  lineage, length(new_names), n_subdiv, spacing))
  
  # ---- add new node coordinates to the shared dp_mst ----
  add_cols <- t(new_coords)                                       # dims x new
  colnames(add_cols) <- new_names
  rownames(add_cols) <- rownames(Y)
  Y_new <- cbind(Y, add_cols)
  cds@principal_graph_aux[[reduction_method]]$dp_mst <- Y_new
  
  # ---- rebuild the densified lineage subgraph ----
  all_nodes <- c(vs, new_names)
  node_df   <- data.frame(name = all_nodes,
                          x = Y_new[1, all_nodes], y = Y_new[2, all_nodes],
                          stringsAsFactors = FALSE)
  cds@graphs[[lineage]] <- graph_from_data_frame(new_edge_df, vertices = node_df,
                                                 directed = FALSE)
  
  # ---- optionally keep the global principal graph consistent ----
  if (isTRUE(update_principal_graph)) {
    G     <- cds@principal_graph[[reduction_method]]
    G_el  <- get.edgelist(G)
    subp  <- el[edge_len > thresh, , drop = FALSE]                # original long edges
    is_sub <- rep(FALSE, nrow(G_el))
    for (i in seq_len(nrow(subp))) {
      a <- subp[i, 1]; b <- subp[i, 2]
      is_sub <- is_sub | (G_el[, 1] == a & G_el[, 2] == b) | (G_el[, 1] == b & G_el[, 2] == a)
    }
    kept        <- G_el[!is_sub, , drop = FALSE]
    chain_edges <- as.matrix(new_edge_df[new_edge_df$from %in% new_names |
                                           new_edge_df$to   %in% new_names, c("from", "to")])
    G_el_new <- rbind(kept, chain_edges)
    all_v    <- union(V(G)$name, new_names)
    vdf      <- data.frame(name = all_v, x = Y_new[1, all_v], y = Y_new[2, all_v],
                           stringsAsFactors = FALSE)
    cds@principal_graph[[reduction_method]] <-
      graph_from_data_frame(as.data.frame(G_el_new, stringsAsFactors = FALSE),
                            vertices = vdf, directed = FALSE)
  }
  
  # ---- re-select cells along the densified graph (fills the gap) ----
  if (isTRUE(reselect_cells)) {
    cds <- isolate_lineage(cds, lineage, sel_clusters = sel_clusters,
                           start_regions = start_regions,
                           starting_clusters = starting_clusters,
                           subset = FALSE, N = N, cl = cl)
  }
  cds
}

cell.selector <- function(path, cells, r, cl){
sel.cells = c()
sel.cells = pbapply(path, 1, selector_sub, cells = cells, r = r, cl = cl, simplify = T)
return(unique(unlist(sel.cells)))
}

cell.selector_sub2 <- function(cell, coords, r){
x2 = cell[1]
y2 = cell[2]
d.x = x2 - coords[1]
d.y = y2 - coords[2]
dist = sqrt(d.x*d.x + d.y*d.y)
if(dist <= r){
return(TRUE)
}
else{
return(FALSE)
}
}

selector_sub <- function(node, cells, r){
x1 = node[1]
y1 = node[2]
res = apply(cells, 1, cell.selector_sub2, coords = c(x1, y1), r = r, simplify = T)
res = names(res[res == TRUE])
return(res)
}

find_start_node <- function(cds){
  nodes = c()
  for(name in names(cds@graphs)){
    sub.graph = cds@graphs[[name]]
    start_end = V(sub.graph)[degree(sub.graph) == 1]$name
    nodes = append(nodes, start_end)
  }
  nodes = as.character(nodes)
  start = names(sort(table(nodes),decreasing=TRUE)[1])
  start
}

