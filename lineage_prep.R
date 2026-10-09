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

find_start_node <- function(cds, which = "subgraph") {
  nodes <- c()
  for (name in names(cds@graphs)) {
    sub.graph <- cds@graphs[[name]]
    if (!inherits(sub.graph, "igraph")) sub.graph <- sub.graph[[which]]   # unwrap converted entries
    start_end <- V(sub.graph)[degree(sub.graph) == 1]$name
    nodes <- append(nodes, start_end)
  }
  nodes <- as.character(nodes)
  names(sort(table(nodes), decreasing = TRUE)[1])
}

isolate_lineage <- function(cds, lineage, sel_clusters = NULL, start_regions = NULL, starting_clusters = NULL,
                            excluded_ages = NULL, excluded_ages_cluster = NULL,
                            subset = FALSE, N = 5, cl = 1, r = 1) {
  stopifnot(is.numeric(r), length(r) == 1, r > 0)
  sel.cells <- isolate_lineage_sub(cds, lineage, sel_clusters = sel_clusters,
                                   start_regions = start_regions, starting_clusters = starting_clusters,
                                   excluded_ages = excluded_ages, excluded_ages_cluster = excluded_ages_cluster,
                                   subset = subset, N = N, cl = cl, r = r)
  cds@lineages[[lineage]] <- sel.cells
  cds
}

isolate_lineage_sub <- function(cds, lineage, sel_clusters = NULL, start_regions = NULL, starting_clusters = NULL,
                                excluded_ages = NULL, excluded_ages_cluster = NULL, subset = FALSE, N = 5, cl = 1, r){
  sub.graph = cds@graphs[[lineage]]
  nodes_UMAP = cds@principal_graph_aux[["UMAP"]]$dp_mst
  if(subset == F){
    nodes_UMAP.sub = as.data.frame(t(nodes_UMAP[,names(V(sub.graph))]))
  }
  else{
    g = principal_graph(cds)[["UMAP"]]
    dd = degree(g)
    names1 = names(dd[dd > 2 | dd == 1])
    names2 = names(dd[dd == 2])
    names2 = sample(names2, length(names2)/subset, replace = F)
    names = c(names1, names2)
    names = intersect(names(V(sub.graph)), names)
    nodes_UMAP.sub = as.data.frame(t(nodes_UMAP[,names]))
  }
  #select cells along the graph
  #mean.dist = path.distance(nodes_UMAP.sub)
  cells_UMAP = as.data.frame(reducedDims(cds)["UMAP"])
  colnames(cells_UMAP) <- toupper(colnames(cells_UMAP))
  cells_UMAP = cells_UMAP[,c("UMAP_1", "UMAP_2")]
  sel.cells = cell.selector(nodes_UMAP.sub, cells_UMAP, r, cl = cl)
  #only keep cells in the progenitor and lineage-specific clusters
  sel.cells1 = c()
  sel.cells2 = sel.cells
  if(length(starting_clusters) > 0){
    sel.cells1 = names(cds@"clusters"[["UMAP"]]$clusters[cds@"clusters"[["UMAP"]]$clusters %in% starting_clusters])
  }
  if(length(start_regions) > 0){
    sel.cells1 = sel.cells1[sel.cells1 %in% rownames(cds@colData[cds@colData$region %in% start_regions,])]
  }
  if(length(sel_clusters) > 0){
    sel.cells2 = names(cds@"clusters"[["UMAP"]]$clusters[cds@"clusters"[["UMAP"]]$clusters %in% sel_clusters])
  }
  
  # Optional: exclude cells in all selected clusters with specific age labels, if no cluster specified, then exclude the cells
  # in all clusters with the specific age labels
  if (length(excluded_ages) > 0) {
    cd <- colData(cds)
    
    # decide which age column the labels belong to (same check as before)
    if (all(excluded_ages %in% cd$age_reorder)) {
      age_col <- "age_reorder"
    } else if (all(excluded_ages %in% cd$age_details)) {
      age_col <- "age_details"
    } else {
      stop("excluded_ages not all found in a single column. Missing from age_reorder: ",
           paste(setdiff(excluded_ages, cd$age_reorder), collapse = ", "),
           " | missing from age_details: ",
           paste(setdiff(excluded_ages, cd$age_details), collapse = ", "))
    }
    
    cluster_labels <- cds@"clusters"[["UMAP"]]$clusters
    
    in_cluster <- if (length(excluded_ages_cluster) > 0) {
      as.character(cluster_labels) %in% as.character(excluded_ages_cluster)
    } else {
      rep(TRUE, length(cluster_labels))     # no cluster given: apply to all cells
    }
    
    cells_to_exclude <- names(cluster_labels)[
      which(in_cluster & cd[names(cluster_labels), age_col] %in% excluded_ages)
    ]
    sel.cells2 <- sel.cells2[!(sel.cells2 %in% cells_to_exclude)]
  }
  cells = unique(c(sel.cells1, sel.cells2))
  sel.cells = sel.cells[sel.cells %in% cells]
  return(sel.cells)
}

get_lineage_object <- function(cds, lineage = FALSE, N = FALSE, recalculate_pt = TRUE){
  start = find_start_node(cds)
  if (lineage != FALSE) {
    sub.graph <- cds@graphs[[lineage]]
    if (is.list(cds@lineages[[lineage]])) {
      sel.cells <- cds@lineages[[lineage]]$name
    } else {
      sel.cells <- cds@lineages[[lineage]]
    }
    if (!is.character(sel.cells)) {
      print("sel cells are not string")
    }
  }
  else{
    sel.cells = colnames(cds)
  }
  sel.cells = sel.cells[sel.cells %in% colnames(cds)]
  nodes_UMAP = cds@principal_graph_aux[["UMAP"]]$dp_mst
  if(N != FALSE){
    if(N < length(sel.cells)){
      sel.cells = sample(sel.cells, N)
    }
  }
  #subset the moncole object
  cds_subset = cds[,sel.cells]
  #set the graph, node and cell UMAP coordinates
  if(lineage == FALSE){
    sub.graph = principal_graph(cds_subset)[["UMAP"]]
  }
  nodes_UMAP <- nodes_UMAP[,names(V(sub.graph))]
  #Reorder the vertices
  degrees <- igraph::degree(sub.graph)
  endpoints <- names(degrees[degrees == 1])
  path_result <- igraph::shortest_paths(sub.graph, from = start, to = endpoints[endpoints != start])
  path_names <- names(path_result$vpath[[1]])   # ordered root -> tip
  # Build the new sequential, zero-padded names
  n <- length(path_names)
  width <- nchar(as.character(n))
  new_names <- paste0("Y_", formatC(seq_len(n), width = width, flag = "0"))
  # Create the old-name -> new-name mapping
  rename_map <- setNames(new_names, path_names)
  igraph::V(sub.graph)$old_name <- igraph::V(sub.graph)$name
  igraph::V(sub.graph)$name <- rename_map[igraph::V(sub.graph)$name]
  #colnames(nodes_UMAP) <- rename_map[colnames(nodes_UMAP)]
  new_colnames <- unname(rename_map[colnames(nodes_UMAP)])
  colnames(nodes_UMAP) <- new_colnames
  cds_subset@principal_graph[["UMAP"]] <- sub.graph
  cds_subset@principal_graph_aux[["UMAP"]]$dp_mst <- nodes_UMAP
  cds_subset@clusters[["UMAP"]]$partitions <- cds_subset@clusters[["UMAP"]]$partitions[colnames(cds_subset)]
  #recalculate closest vertex and pseudotime for the selected cells
  if(recalculate_pt == TRUE){
    source_url("https://raw.githubusercontent.com/cole-trapnell-lab/monocle3/master/R/learn_graph.R")
    cds_subset <- project2MST(cds_subset, project_point_to_line_segment, F, T, "UMAP", nodes_UMAP)
    cds_subset <- order_cells(cds_subset, root_pr_nodes = unname(rename_map[start]))
  }
  return(cds_subset)
}

pt_recalculate <- function(cds, lineages = names(cds@lineages), N = FALSE, recalculate_pt = TRUE) {
  for (lineage in lineages) {
    message("Processing lineage ", lineage)
    # skip lineages already converted by an earlier run
    sub.graph <- cds@graphs[[lineage]]
    if (!inherits(sub.graph, "igraph")) {
      message("  skipping ", lineage, ": cds@graphs[[lineage]] is not an igraph (already converted?)")
      next
    }
    lin <- cds@lineages[[lineage]]
    cell_names <- if (is.list(lin)) lin$name else lin
    lineage_sub <- get_lineage_object(cds, lineage = lineage, N = N, recalculate_pt = recalculate_pt)
    new.graph <- lineage_sub@principal_graph$UMAP
    dp_mst    <- lineage_sub@principal_graph_aux@listData[["UMAP"]][["dp_mst"]]
    
    # consistency checks, done before anything is stored
    n_sub <- igraph::vcount(sub.graph)
    n_new <- igraph::vcount(new.graph)
    n_dp  <- ncol(dp_mst)
    if (!(n_sub == n_new && n_new == n_dp)) {
      stop("lineage '", lineage, "': sizes do not match -> subgraph: ", n_sub,
           " vertices, reordered graph: ", n_new, " vertices, dp_mst: ", n_dp, " columns", call. = FALSE)
    }
    if (!setequal(igraph::V(new.graph)$name, colnames(dp_mst))) {
      stop("lineage '", lineage, "': vertex names of the reordered graph and the dp_mst columns differ",
           call. = FALSE)
    }
    
    pt <- pseudotime(lineage_sub)
    stopifnot(all(cell_names %in% names(pt)))
    
    cds@lineages[[lineage]] <- list(
      name = cell_names,
      updated_pt = pt[cell_names]
    )
    cds@graphs[[lineage]] <- list(
      subgraph = sub.graph,
      subgraph_reorder = new.graph,
      dp_mst = dp_mst
    )
    rm(lineage_sub, pt, sub.graph, new.graph, dp_mst); gc()
  }
  cds
}
