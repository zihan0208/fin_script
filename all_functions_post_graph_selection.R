densify_all_lineages <- function(cds, lineages = names(cds@graphs), reselect_cells = FALSE, ...) {
  missing <- setdiff(lineages, names(cds@graphs))
  if (length(missing) > 0) stop("lineages not found in cds@graphs: ", paste(missing, collapse = ", "), call. = FALSE)
  
  for (lineage in lineages) {
    cds <- densify_lineage_graph(cds, lineage, reselect_cells = reselect_cells, ...)
  }
  cds
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


isolate_lineage <- function(cds, lineage, sel_clusters = NULL, start_regions = NULL, starting_clusters = NULL,
                            excluded_ages = NULL, excluded_ages_cluster = NULL,
                            subset = FALSE, cl = 1, r = 1) {
  stopifnot(is.numeric(r), length(r) == 1, r > 0)
  sel.cells <- isolate_lineage_sub(cds, lineage, sel_clusters = sel_clusters,
                                   start_regions = start_regions, starting_clusters = starting_clusters,
                                   excluded_ages = excluded_ages, excluded_ages_cluster = excluded_ages_cluster,
                                   subset = subset, cl = cl, r = r)
  cds@lineages[[lineage]] <- sel.cells
  cds
}


isolate_lineage_sub <- function(cds, lineage, sel_clusters = NULL, start_regions = NULL, starting_clusters = NULL,
                                excluded_ages = NULL, excluded_ages_cluster = NULL, subset = FALSE, cl = 1, r){
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


compress_lineage <- function(cds, lineage, n_metacells = 1000,
                             pt_name = "updated_pt") {
  
  lin <- cds@lineages[[lineage]]
  if (!is.list(lin) || is.null(lin$name) || is.null(lin[[pt_name]]))
    stop("lineage '", lineage, "' needs a list with $name and $", pt_name)
  
  pt_vec <- lin[[pt_name]]
  if (is.null(names(pt_vec))) {
    if (length(pt_vec) != length(lin$name)) stop("unnamed pseudotime does not match $name in length")
    names(pt_vec) <- lin$name
  }
  
  # cells that exist in the object and have a pseudotime
  cells <- lin$name[lin$name %in% colnames(cds)]
  cells <- cells[cells %in% names(pt_vec) & !is.na(pt_vec[cells])]
  if (length(cells) < n_metacells) warning(lineage, ": fewer cells than metacells (", length(cells), ")")
  
  cds_subset = cds[,cells]
  #preprare raw count matrix
  exp = as.data.frame(as.matrix(exprs(cds_subset)))
  exp_sum <- t(exp)
  exp_sum = exp_sum[,rownames(cds)]
  #prepare size factor
  size_factor <- (pData(cds_subset)[, 'Size_Factor'])
  size_factor <- size_factor[rownames(exp_sum)]
  #prepare pseudotime
  #pt <- cds@lineages[[lineage]][['updated_pt']]
  #pt <- cds@lineages[[lineage]]$gap_reduced_pt
  pt <- as.data.frame(pt_vec)
  pt <- pt[rownames(exp_sum), , drop = FALSE]
  colnames(pt) <- c("pseudotime")
  #prepare umap
  UMAP <- reducedDims(cds_subset)[["UMAP"]]
  UMAP <- UMAP[rownames(exp_sum),]
  exp_sum = cbind(pt, UMAP, size_factor, exp_sum)
  exp_sum = exp_sum[order(exp_sum$pseudotime),]
  exp_sum$meta_cell <- cut(rank(exp_sum$pseudotime), breaks = n_metacells, labels = FALSE)
  gene_cols <- setdiff(
    colnames(exp_sum),
    c("pseudotime", "umap_1", "umap_2", "size_factor", "meta_cell")
  )
  print(paste0("Compressing lineage ", lineage, " with sum"))
  meta_sum <- exp_sum %>%
    group_by(meta_cell) %>%
    summarise(
      n_cells = n(),                     # number of cells in this meta-cell
      pseudotime = mean(pseudotime),     # mean pseudotime
      umap_1 = mean(umap_1),             # mean UMAP_x
      umap_2 = mean(umap_2),             # mean UMAP_y
      size_factor = sum(size_factor),    # sum size factors
      across(all_of(gene_cols), sum)   # sum gene counts
    ) %>%
    ungroup()
  meta_sum_ordered <- meta_sum[order(meta_sum$pseudotime), ]
  
  exp_mean <- exp
  exp_mean = (t(exp_mean)) /  (pData(cds_subset)[, 'Size_Factor'])
  exp_mean = exp_mean[,rownames(cds)]
  pt <- pt[rownames(exp_mean), , drop = FALSE]
  UMAP <- UMAP[rownames(exp_mean),]
  exp_mean = cbind(pt, UMAP, exp_mean)
  exp_mean = exp_mean[order(exp_mean$pseudotime),]
  exp_mean$meta_cell <- cut(rank(exp_mean$pseudotime), breaks = n_metacells, labels = FALSE)
  gene_cols <- setdiff(
    colnames(exp_mean),
    c("pseudotime", "umap_1", "umap_2", "meta_cell")
  )
  print(paste0("Compressing lineage ", lineage, " with mean"))
  meta_mean <- exp_mean %>%
    group_by(meta_cell) %>%
    summarise(
      n_cells = n(),                     # number of cells in this meta-cell
      pseudotime = mean(pseudotime),     # mean pseudotime
      umap_1 = mean(umap_1),             # mean UMAP_x
      umap_2 = mean(umap_2),             # mean UMAP_y
      across(all_of(gene_cols), mean)   # mean gene counts
    ) %>%
    ungroup()
  meta_mean_ordered <- meta_mean[order(meta_mean$pseudotime), ]
  list(sum = meta_sum_ordered, mean = meta_mean_ordered)
}

compress_all_lineages <- function(cds, lineages = names(cds@lineages), save_dir = NULL, ...) {
  if (!is.null(save_dir)) dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)
  
  for (lineage in lineages) {
    message("Compressing lineage ", lineage)
    res <- tryCatch(compress_lineage(cds, lineage, ...),
                    error = function(e) { warning("skipping ", lineage, ": ", conditionMessage(e)); NULL })
    if (is.null(res)) next
    
    if (is.null(cds@expression[[lineage]])) cds@expression[[lineage]] <- list()
    cds@expression[[lineage]]$sum  <- res$sum
    cds@expression[[lineage]]$mean <- res$mean
    
    if (!is.null(save_dir))
      saveRDS(res, file.path(save_dir, paste0(gsub("[^A-Za-z0-9_.-]", "_", lineage), "_meta.rds")))
    rm(res); gc()
  }
  cds
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

plot_combined_graph <- function(cds, reduction_method = "UMAP", color_cells_by = NULL, alpha = 0.4,
                                label_clusters = FALSE, highlight_clusters = NULL,
                                other_color = "grey85", show_others = TRUE,
                                highlight_lineage = NULL, highlight_color = "orange",
                                highlight_size = 0.3, highlight_alpha = 0.6,
                                highlight_on_top = TRUE, graph_which = "subgraph",
                                highlight_vertices = NULL, highlight_cells = NULL,
                                point_color = "red", point_size = 3, label_points = TRUE) {
  library(igraph)
  library(ggplot2)
  library(dplyr)
  
  # accept both forms: a plain igraph, or the converted list(subgraph, subgraph_reorder, dp_mst)
  graph_list <- lapply(cds@graphs, function(g) {
    if (inherits(g, "igraph")) g else g[[graph_which]]
  })
  stopifnot(all(vapply(graph_list, inherits, logical(1), what = "igraph")))
  
  combined_graph <- Reduce(igraph::union, graph_list)
  graph_nodes <- V(combined_graph)$name
  
  node_coords <- as.data.frame(t(cds@principal_graph_aux[[reduction_method]]$dp_mst))
  node_coords$node <- rownames(node_coords)
  colnames(node_coords)[1:2] <- c("x", "y")
  node_coords <- node_coords %>% filter(node %in% graph_nodes)
  
  el <- igraph::get.edgelist(combined_graph)
  edges_df <- data.frame(
    from = el[,1], to = el[,2],
    x    = node_coords$x[match(el[,1], node_coords$node)],
    y    = node_coords$y[match(el[,1], node_coords$node)],
    xend = node_coords$x[match(el[,2], node_coords$node)],
    yend = node_coords$y[match(el[,2], node_coords$node)]
  )
  
  cell_coords <- as.data.frame(reducedDims(cds)[[reduction_method]], stringsAsFactors = FALSE)
  colnames(cell_coords) <- c("x", "y")
  cell_coords$cell_id <- rownames(cell_coords)
  has_color <- !is.null(color_cells_by) && color_cells_by %in% colnames(colData(cds))
  if (has_color) {
    cell_coords$color <- as.character(colData(cds)[cell_coords$cell_id, color_cells_by])
  } else {
    cell_coords$color <- NA_character_
  }
  
  if (!is.null(highlight_clusters) && !has_color)
    stop("highlight_clusters needs color_cells_by (the column that holds the cluster labels)")
  
  highlight_coords <- NULL
  if (!is.null(highlight_lineage)) {
    if (!highlight_lineage %in% names(cds@lineages))
      stop("lineage '", highlight_lineage, "' not found in cds@lineages")
    lin <- cds@lineages[[highlight_lineage]]
    lin_cells <- if (is.list(lin)) lin$name else lin
    highlight_coords <- cell_coords[cell_coords$cell_id %in% lin_cells, ]
    if (nrow(highlight_coords) == 0) stop("none of the lineage cells were found in the UMAP")
  }
  
  hl_layer <- NULL
  if (!is.null(highlight_coords)) {
    hl_layer <- geom_point(data = highlight_coords, aes(x = x, y = y),
                           color = highlight_color, size = highlight_size, alpha = highlight_alpha)
  }
  
  vertex_pts <- NULL
  if (!is.null(highlight_vertices)) {
    dp   <- cds@principal_graph_aux[[reduction_method]]$dp_mst
    miss <- setdiff(highlight_vertices, colnames(dp))
    if (length(miss) > 0) warning(length(miss), " vertices not found in dp_mst: ", paste(head(miss, 5), collapse = ", "))
    keep <- intersect(highlight_vertices, colnames(dp))
    vertex_pts <- data.frame(name = keep, x = dp[1, keep], y = dp[2, keep])
  }
  
  cell_pts <- NULL
  if (!is.null(highlight_cells)) {
    cell_pts <- cell_coords[cell_coords$cell_id %in% highlight_cells, ]
    if (nrow(cell_pts) < length(highlight_cells))
      warning(length(highlight_cells) - nrow(cell_pts), " cells not found in the UMAP")
  }
  
  p <- ggplot()
  if (!highlight_on_top) p <- p + hl_layer
  
  if (!has_color) {
    # plain grey cells, no color aesthetic, so there is no "NA" legend
    p <- p + geom_point(data = cell_coords, aes(x = x, y = y),
                        color = other_color, size = 0.01, alpha = alpha)
    label_df <- NULL
  } else if (is.null(highlight_clusters)) {
    p <- p + geom_point(data = cell_coords, aes(x = x, y = y, color = color), size = 0.01, alpha = alpha)
    label_df <- cell_coords
  } else {
    is_hl <- cell_coords$color %in% as.character(highlight_clusters)
    if (!any(is_hl)) stop("none of highlight_clusters found in ", color_cells_by)
    if (show_others) {
      p <- p + geom_point(data = cell_coords[!is_hl, ], aes(x = x, y = y),
                          color = other_color, size = 0.01, alpha = alpha)
    }
    p <- p + geom_point(data = cell_coords[is_hl, ], aes(x = x, y = y, color = color),
                        size = 0.01, alpha = alpha)
    label_df <- cell_coords[is_hl, ]
  }
  
  if (highlight_on_top) p <- p + hl_layer          # above the cells, below the trajectory
  
  p <- p +
    geom_segment(data = edges_df, aes(x = x, y = y, xend = xend, yend = yend), color = "black", size = 1) +
    geom_point(data = node_coords, aes(x = x, y = y), color = "blue", size = 0.2) +
    theme_minimal() +
    coord_fixed()
  
  if (!is.null(vertex_pts) && nrow(vertex_pts) > 0) {
    p <- p + geom_point(data = vertex_pts, aes(x = x, y = y), shape = 21, fill = "yellow",
                        color = point_color, size = point_size, stroke = 1.2)
    if (label_points)
      p <- p + ggrepel::geom_text_repel(data = vertex_pts, aes(x = x, y = y, label = name),
                                        size = 3.5, fontface = "bold", bg.color = "white", bg.r = 0.15, max.overlaps = Inf)
  }
  if (!is.null(cell_pts) && nrow(cell_pts) > 0) {
    p <- p + geom_point(data = cell_pts, aes(x = x, y = y), color = point_color, size = point_size - 1)
  }
  
  if (label_clusters && !is.null(label_df)) {
    centroids <- label_df %>%
      group_by(color) %>%
      summarise(x = median(x), y = median(y), .groups = "drop")
    p <- p + ggrepel::geom_text_repel(data = centroids, aes(x = x, y = y, label = color),
                                      size = 5, color = "black", fontface = "bold",
                                      bg.color = "white", bg.r = 0.15, max.overlaps = Inf)
  }
  
  p
}

get_branch_points <- function(graph_list, full_graph, graph_which = "subgraph") {
  full_names <- igraph::V(full_graph)$name
  # one igraph per lineage, whichever form is stored
  graphs <- lapply(graph_list, function(g) if (inherits(g, "igraph")) g else g[[graph_which]])
  stopifnot(all(vapply(graphs, inherits, logical(1), what = "igraph")))
  combined_graph <- NULL
  for (g in graphs) {
    keep <- intersect(igraph::V(g)$name, full_names)
    g_orig <- igraph::induced_subgraph(full_graph, igraph::V(full_graph)[name %in% keep])
    combined_graph <- if (is.null(combined_graph)) g_orig else igraph::union(combined_graph, g_orig)
  }
  combined_graph <- igraph::simplify(combined_graph)
  branch_points <- igraph::V(combined_graph)[igraph::degree(combined_graph) >= 3]
  branch_points$name
}

align_branch_points_lineage <- function(cds, lineage, branch_points,
                                        shift_to_zero = TRUE, umap_name = "UMAP") {
  
  # graph with renamed vertices, plus its own coordinate matrix
  entry <- cds@graphs[[lineage]]
  if (inherits(entry, "igraph") || is.null(entry$subgraph_reorder) || is.null(entry$dp_mst))
    stop("cds@graphs[['", lineage, "']] must be a list with $subgraph_reorder and $dp_mst")
  
  g <- entry$subgraph_reorder
  if (is.null(igraph::V(g)$old_name)) stop("subgraph_reorder has no 'old_name' vertex attribute")
  if (anyDuplicated(igraph::V(g)$old_name)) stop("old_name is not unique, the name map would be ambiguous")
  
  # original branch point name -> this lineage's new vertex name (NA if the lineage lacks it)
  name_map  <- setNames(igraph::V(g)$name, igraph::V(g)$old_name)
  bp_new    <- unname(name_map[branch_points])
  
  # pseudotime of the lineage's cells
  lin <- cds@lineages[[lineage]]
  if (!is.list(lin) || is.null(lin$name) || is.null(lin$updated_pt))
    stop("cds@lineages[['", lineage, "']] must be a list with $name and $updated_pt")
  pt <- lin$updated_pt
  if (is.null(names(pt))) names(pt) <- lin$name
  pt <- pt[!is.na(pt)]
  if (shift_to_zero) pt <- pt - min(pt)
  
  cells <- intersect(lin$name, intersect(names(pt), colnames(cds)))
  if (length(cells) == 0) stop("no lineage cells with a pseudotime were found in the object")
  
  # coordinates, taken without subsetting the whole cds
  cell_umap   <- as.matrix(reducedDims(cds)[[umap_name]][cells, 1:2, drop = FALSE])
  vertex_umap <- t(entry$dp_mst)
  
  closest <- setNames(rep(NA_character_, length(branch_points)), branch_points)
  target  <- setNames(rep(NA_real_,      length(branch_points)), branch_points)
  
  for (i in seq_along(branch_points)) {
    v <- bp_new[i]
    if (is.na(v) || !(v %in% rownames(vertex_umap))) next        # branch point not in this lineage
    v_coord <- vertex_umap[v, 1:2]
    d  <- sqrt((cell_umap[, 1] - v_coord[1])^2 + (cell_umap[, 2] - v_coord[2])^2)
    cc <- cells[which.min(d)]
    closest[i] <- cc
    target[i]  <- unname(pt[cc])
  }
  
  res <- rbind(closest_cells = closest, target_pts = target)
  colnames(res) <- branch_points
  res
}
align_branch_points <- function(cds, branch_points, lineages = names(cds@lineages), ...) {
  results_list <- list()
  skipped <- character(0)
  
  for (lineage in lineages) {
    message("Aligning lineage ", lineage)
    res <- tryCatch(align_branch_points_lineage(cds, lineage, branch_points, ...),
                    error = function(e) { warning("skipping ", lineage, ": ", conditionMessage(e)); NULL })
    if (is.null(res)) { skipped <- c(skipped, lineage); next }
    results_list[[lineage]] <- res
  }
  
  if (length(skipped) > 0) message("Skipped: ", paste(skipped, collapse = ", "))
  results_list
}

collect_by_branch_point <- function(results_list, branch_points) {
  stopifnot(length(results_list) > 0)
  
  final_output <- list()
  for (bp in branch_points) {
    bp_data <- sapply(results_list, function(x) {
      if (bp %in% colnames(x)) x[, bp] else c(NA_character_, NA_character_)   # branch point missing from this lineage
    })
    rownames(bp_data) <- c("closest_cell", "target_pseudotime")
    final_output[[paste0("BP_", bp)]] <- bp_data
  }
  final_output
}

get_fixed_targets <- function(cds, final_output, expr_slot = "sum", pt_col = "pseudotime", target_row = 2) {
  
  stopifnot(length(final_output) > 0, !is.null(names(final_output)))
  
  # pseudotime column of a lineage's metacell table: by name, falling back to position 3
  pt_of <- function(lin) {
    tab <- cds@expression[[lin]][[expr_slot]]
    if (is.null(tab)) return(NULL)
    col <- if (pt_col %in% colnames(tab)) pt_col else 3
    as.numeric(unlist(tab[, col]))
  }
  
  fixed_targets <- setNames(rep(NA_real_, length(final_output)), names(final_output))
  
  for (bp in names(final_output)) {
    m <- final_output[[bp]]
    raw_pts <- as.numeric(m[target_row, ])
    names(raw_pts) <- colnames(m)
    active_lins <- names(raw_pts)[!is.na(raw_pts)]
    if (length(active_lins) == 0) next
    
    # lineages that have an expression table
    pts <- lapply(active_lins, pt_of); names(pts) <- active_lins
    pts <- pts[!vapply(pts, is.null, logical(1))]
    if (length(pts) == 0) next
    
    # 1. span of every lineage that shares this breakpoint
    spans <- vapply(pts, function(x) diff(range(x, na.rm = TRUE)), numeric(1))
    
    # 2. the "master" (longest) lineage for this breakpoint
    master <- names(spans)[which.max(spans)]
    max_span <- spans[[master]]
    if (!is.finite(max_span) || max_span <= 0) next
    
    # 3. position of the breakpoint on the master's 0-1 axis
    start_master <- min(pts[[master]], na.rm = TRUE)
    fixed_targets[bp] <- (raw_pts[[master]] - start_master) / max_span
  }
  
  sort(fixed_targets[!is.na(fixed_targets)])
}

store_aligned_pseudotime <- function(cds, bp_aligned_list, slot_prefix = "aligned_",
                                     store = "streched", overwrite = TRUE) {
  
  stopifnot(is.list(bp_aligned_list), length(bp_aligned_list) > 0, !is.null(names(bp_aligned_list)))
  
  for (bp_name in names(bp_aligned_list)) {
    per_lineage <- bp_aligned_list[[bp_name]]
    for (lin in names(per_lineage)) {
      slot_name <- paste0(slot_prefix, bp_name)
      
      # make sure the container exists
      if (is.null(cds@pseudotime[[lin]])) cds@pseudotime[[lin]] <- list()
      if (is.null(cds@pseudotime[[lin]][[store]]) || !is.list(cds@pseudotime[[lin]][[store]]))
        cds@pseudotime[[lin]][[store]] <- list()
      
      if (!overwrite && !is.null(cds@pseudotime[[lin]][[store]][[slot_name]])) next
      cds@pseudotime[[lin]][[store]][[slot_name]] <- per_lineage[[lin]]
    }
  }
  cds
}


get_bp_specific_alignments <- function(cds, lineages, bp_results, fixed_targets) {
  # This list will store: results[["BP_NAME"]][["LINEAGE_NAME"]]
  all_alignments <- list()
  # Pre-calculate lineage boundaries
  lin_info <- lapply(lineages, function(lin) {
    vals <- as.numeric(unlist(cds@expression[[lin]]$sum[, 3]))
    list(min = min(vals), max = max(vals), span = diff(range(vals)))
  })
  names(lin_info) <- lineages
  
  # Iterate through each Branch Point
  for (bp_name in names(fixed_targets)) {
    target_pos <- fixed_targets[bp_name]
    all_alignments[[bp_name]] <- list()
    # Identify which lineages have this BP
    # Accessing the second row (target_pseudotime) for all lineages
    raw_pts_for_bp <- bp_results[[bp_name]][2, ]
    valid_lins <- names(raw_pts_for_bp)[!is.na(raw_pts_for_bp)]
    
    for (lin in valid_lins) {
      # Extract raw meta-cell PT
      pt_meta <- as.numeric(unlist(cds@expression[[lin]]$sum[, 3]))
      names(pt_meta) <- rownames(cds@expression[[lin]]$sum)
      
      # The raw anchor for this lineage
      raw_anchor <- as.numeric(raw_pts_for_bp[lin])
      
      # Define 2 segments: [Start to BP] and [BP to End]
      raw_bounds <- c(lin_info[[lin]]$min, raw_anchor, lin_info[[lin]]$max)
      new_bounds <- c(0, target_pos, 1.0)
      
      pt_norm <- pt_meta
      for (i in 1:2) {
        lr <- raw_bounds[i]; ur <- raw_bounds[i+1]
        ln <- new_bounds[i]; un <- new_bounds[i+1]
        
        # Segment mask logic [lower, upper) for first, [lower, upper] for last
        mask <- if(i == 1) (pt_meta >= lr & pt_meta < ur) else (pt_meta >= lr & pt_meta <= ur)
        
        if (any(mask) && (ur > lr)) {
          pt_norm[mask] <- (pt_meta[mask] - lr) / (ur - lr) * (un - ln) + ln
        }
      }
      # Store the result for this specific lineage under this specific BP
      all_alignments[[bp_name]][[lin]] <- pt_norm
    }
  }
  return(all_alignments)
}



plot_lineage_alignment <- function(cds, lineages, gene, bp_id = NULL, 
                                   overlay = FALSE,
                                   targets = fixed_targets,
                                   custom_lineage_colors = NULL) {
  # 1. Collect data
  plot_df <- do.call(rbind, lapply(lineages, function(lin) {
    get_plotting_data(cds, lin, gene, bp_id)
  }))
  
  if (is.null(plot_df) || nrow(plot_df) == 0) stop("No data found.")
  
  # --- Lineage Ordering (Numerical) ---
  all_lineage_names <- names(cds@lineages)
  nums <- as.numeric(gsub("[^0-9]", "", all_lineage_names))
  ordered_lineages <- if (!all(is.na(nums))) all_lineage_names[order(nums)] else sort(all_lineage_names)
  
  # Colors
  n_lin <- length(ordered_lineages)
  base_colors <- brewer.pal(9, "Blues")[-c(1,2,3)]
  lineage_cols <- colorRampPalette(base_colors)(n_lin)
  names(lineage_cols) <- ordered_lineages
  
  if (!is.null(custom_lineage_colors)) {
    for (lin in names(custom_lineage_colors)) {
      if (lin %in% names(lineage_cols)) lineage_cols[lin] <- custom_lineage_colors[lin]
    }
  }
  lineage_cols_lines <- darken(lineage_cols, amount = 0.3)
  lineage_cols_points <- rev(lineage_cols)
  names(lineage_cols_points) <- ordered_lineages
  
  # Set factor levels for legend order
  plot_df$Lineage <- factor(plot_df$Lineage, levels = ordered_lineages)
  
  # 2. Numeric Sorting for Axis Factor
  unique_axes <- unique(as.character(plot_df$Axis))
  ax_nums <- as.numeric(gsub("[^0-9]", "", unique_axes))
  if (!all(is.na(ax_nums))) {
    plot_df$Axis <- factor(plot_df$Axis, levels = unique_axes[order(ax_nums)])
  }
  
  # 3. Construct Dynamic Title
  clean_bp_names <- if(!is.null(bp_id)) sub("aligned_", "", bp_id) else "Raw"
  if (is.null(bp_id)) {
    title_text <- paste("Alignment:", gene, "| Raw Pseudotime")
  } else if (identical(bp_id, "all")) {
    title_text <- paste("Gene:", gene, "| All Branch Points")
  } else {
    title_text <- paste0("Gene: ", gene, " | BP: ", paste(clean_bp_names, collapse = ", "))
  }
  
  # 4. Prepare Vertical Line Data
  v_lines <- data.frame(
    Axis = sub("aligned_", "", names(targets)),
    TargetPos = as.numeric(targets)
  )
  
  # 5. Build Base Plot
  df_points <- subset(plot_df, Type == "Raw")
  df_points$color_mirror <- lineage_cols_points[df_points$Lineage]
  
  p <- ggplot(plot_df, aes(x = PT, y = Expression)) +
    geom_point(data = df_points, 
               aes(color = Lineage),
               color = df_points$color_mirror,
               size = 0.8, alpha = 0.3) +
    geom_line(data = subset(plot_df, Type == "Fit"), 
              aes(color = Lineage),
              linewidth = 1.2) +
    scale_color_manual(values = lineage_cols) +
    labs(title = title_text,
         x = "Stretched Pseudotime (0-1)",
         y = "log(Expression + 1)") +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "right",
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
    )
  
  # 6. Refined Faceting & Overlay Logic
  if (!is.null(bp_id)) {
    v_lines_sub <- v_lines[v_lines$Axis %in% levels(plot_df$Axis), ]
    num_axes <- length(levels(plot_df$Axis))
    
    if (!overlay) {
      if (num_axes > 1) {
        # Scenario 1: Multiple BPs selected, show them in separate BP boxes
        p <- p + facet_wrap(~Axis)
      } else {
        # Scenario 2: Single BP selected, show each Lineage in its own box
        p <- p + facet_wrap(~Lineage)
      }
      
      # Add vertical lines to the facets
      p <- p + geom_vline(data = v_lines_sub, aes(xintercept = TargetPos), 
                          linetype = "dashed", color = "black", alpha = 0.5)
      
    } else {
      # Scenario 3: Overlay is TRUE (all BPs and all Lineages on one plot)
      p <- p + geom_vline(data = v_lines_sub, aes(xintercept = TargetPos), 
                          linetype = "dashed", color = "black", alpha = 0.4)
      
      if (num_axes > 1) {
        p <- p + labs(subtitle = "Overlay Mode: Multiple Branch Points")
      }
    }
  }
  
  
  return(p)
}

plot_lineage_alignment_raw <- function(cds, lineages, gene, bp_id = NULL, 
                                       overlay = FALSE,
                                       targets = fixed_targets,
                                       custom_lineage_colors = NULL) {
  # 1. Collect data (Assuming get_plotting_data still returns the 'Raw' rows)
  plot_df <- do.call(rbind, lapply(lineages, function(lin) {
    get_plotting_data(cds, lin, gene, bp_id)
  }))
  
  if (is.null(plot_df) || nrow(plot_df) == 0) stop("No data found.")
  
  # Filter to ensure we only have Raw points (in case get_plotting_data still returns Fit)
  plot_df <- subset(plot_df, Type == "Raw")
  
  # --- Numerical Lineage Ordering for Legend ---
  present_lineages <- unique(as.character(plot_df$Lineage))
  nums <- as.numeric(gsub("[^0-9]", "", present_lineages))
  ordered_lineages <- if (!all(is.na(nums))) present_lineages[order(nums)] else sort(present_lineages)
  
  # Colors
  n_lin <- length(ordered_lineages)
  base_colors <- brewer.pal(9, "Blues")[-c(1,2,3)]
  lineage_cols <- colorRampPalette(base_colors)(n_lin)
  names(lineage_cols) <- ordered_lineages
  
  if (!is.null(custom_lineage_colors)) {
    for (lin in names(custom_lineage_colors)) {
      if (lin %in% names(lineage_cols)) lineage_cols[lin] <- custom_lineage_colors[lin]
    }
  }
  
  # Set factor levels for numerical sorting in legend
  plot_df$Lineage <- factor(plot_df$Lineage, levels = ordered_lineages)
  
  # 2. Numerical Sorting for Axis Factor (Branch Points)
  unique_axes <- unique(as.character(plot_df$Axis))
  ax_nums <- as.numeric(gsub("[^0-9]", "", unique_axes))
  if (!all(is.na(ax_nums))) {
    plot_df$Axis <- factor(plot_df$Axis, levels = unique_axes[order(ax_nums)])
  }
  
  # 3. Dynamic Title
  clean_bp_names <- if(!is.null(bp_id)) sub("aligned_", "", bp_id) else "Raw"
  title_text <- if (identical(bp_id, "all")) {
    paste("Raw Alignment:", gene, "| All Branches")
  } else {
    paste0("Raw Alignment: ", gene, " | BP: ", paste(clean_bp_names, collapse = ", "))
  }
  
  # 4. Vertical Line Data
  v_lines <- data.frame(
    Axis = sub("aligned_", "", names(targets)),
    TargetPos = as.numeric(targets)
  )
  
  # 5. Build Plot (Points Only)
  p <- ggplot(plot_df, aes(x = PT, y = Expression, color = Lineage)) +
    geom_point(size = 0.8, alpha = 0.3) +
    scale_color_manual(values = lineage_cols) +
    labs(title = title_text,
         x = "Stretched Pseudotime (0-1)",
         y = "log(Expression + 1)") +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "right",
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
    )
  
  # 6. Faceting Logic
  if (!is.null(bp_id)) {
    v_lines_sub <- v_lines[v_lines$Axis %in% levels(plot_df$Axis), ]
    num_axes <- length(levels(plot_df$Axis))
    
    if (!overlay) {
      # Multi-branch: facet by Branch. Single-branch: facet by Lineage.
      if (num_axes > 1) { p <- p + facet_wrap(~Axis) } 
      else { p <- p + facet_wrap(~Lineage) }
      
      p <- p + geom_vline(data = v_lines_sub, aes(xintercept = TargetPos), 
                          linetype = "dashed", color = "black", alpha = 0.5)
    } else {
      # Overlay everything in one box
      p <- p + geom_vline(data = v_lines_sub, aes(xintercept = TargetPos), 
                          linetype = "dashed", color = "black", alpha = 0.4)
    }
  }
  
  return(p)
}

get_plotting_data <- function(cds, lineage, gene, bp_id = NULL) {
  # 1. Get Raw Expression (Logged)
  expr_val <- log(cds@expression[[lineage]]$mean[[gene]] + 1)
  
  # Helper to create the 1000-point fit data frame
  get_fit_df <- function(lin, bp_slot) {
    # Check if the lineage and branch point slot exist in expectation
    # 1. Is the lineage even in the expectation list?
    if (!(lin %in% names(cds@expectation))) return(NULL)
    
    # 2. Get the branch names for this lineage
    bp_names <- names(cds@expectation[[lin]])
    
    # 3. If bp_names is NULL or the specific slot isn't there, STOP and return NULL
    if (is.null(bp_names) || !(bp_slot %in% bp_names)) {
      return(NULL)
    }
    
    # 4. Now it is safe to check the gene
    gene_data <- cds@expectation[[lin]][[bp_slot]]
    if (!(gene %in% names(gene_data))) return(NULL)
    
    fit <- gene_data[[gene]]$prediction
    if (is.null(fit)) return(NULL)
    
    clean_bp <- sub("aligned_", "", bp_slot)
    return(data.frame(
      PT = seq(0, 1, length.out = 1000),
      Expression = log(fit + 1),
      Lineage = lin,
      Axis = clean_bp,
      Type = "Fit"
    ))
  }
  
  # 2. Case A: No BP specified (Raw Scaled PT)
  if (is.null(bp_id)) {
    pt_val <- as.numeric(cds@pseudotime[[lineage]]$scaled)
    return(data.frame(
      PT = pt_val, 
      Expression = expr_val, 
      Lineage = lineage, 
      Axis = "Evenly scaled Pseudotime",
      Type = "Raw"
    ))
  }
  
  # 3. Case B: "all" branch points
  if (bp_id == "all") {
    all_slots <- names(cds@pseudotime[[lineage]]$streched)
    if (length(all_slots) == 0) return(NULL)
    
    combined_list <- lapply(all_slots, function(slot) {
      clean_bp <- sub("aligned_", "", slot)
      df_raw <- data.frame(
        PT = cds@pseudotime[[lineage]]$streched[[slot]],
        Expression = expr_val,
        Lineage = lineage,
        Axis = clean_bp,
        Type = "Raw"
      )
      # If fit is NULL, only df_raw is returned via rbind
      df_fit <- get_fit_df(lineage, slot)
      return(rbind(df_raw, df_fit))
    })
    return(do.call(rbind, combined_list))
    
  } else {
    # 4. Case C: One specific BP
    slot_name <- if(grepl("aligned_", bp_id)) bp_id else paste0("aligned_", bp_id)
    pt_val <- cds@pseudotime[[lineage]]$streched[[slot_name]]
    clean_bp <- sub("aligned_", "", slot_name)
    
    if (is.null(pt_val)) return(NULL)
    
    df_raw <- data.frame(
      PT = pt_val, 
      Expression = expr_val, 
      Lineage = lineage, 
      Axis = clean_bp, 
      Type = "Raw"
    )
    
    # If get_fit_df returns NULL, only df_raw is kept
    df_fit <- get_fit_df(lineage, slot_name)
    
    return(rbind(df_raw, df_fit))
  }
}

calculate_separate_wald_end <- function(fit_A, fit_B, l2fc = 0.4, eigenThresh = 1e-8, interior_knots, n_knots) {
  # 1. Extract betas and vcovs
  betaA <- fit_A$beta
  betaB <- fit_B$beta
  VA <- fit_A$vcov
  VB <- fit_B$vcov
  L <- fit_A$X_end_2                     # <-- swapped
  
  # 2. Calculate the difference (Contrast)
  estFC <- L %*% betaA - L %*% betaB
  logFCCutoff <- log(2^l2fc)
  est <- sign(estFC) * pmax(0, abs(estFC) - logFCCutoff)
  
  # 3. Calculate Combined Sigma (Independence assumption)
  sigma <- L %*% (VA + VB) %*% t(L)
  
  # 4. Eigen-decomposition
  eSigma <- eigen(sigma, symmetric = TRUE)
  
  # 5. Determine rank
  r <- try(sum(eSigma$values / eSigma$values[1] > eigenThresh), silent = TRUE)
  
  if (inherits(r, "try-error") || is.na(r) || r < 1) {
    return(c(waldStat = NA, df = NA, p_val = NA, log2FC = NA, predictA = NA, predictB = NA))
  }
  
  if (r == 1) {
    halfCovInv <- eSigma$vectors[, 1, drop = FALSE] * (1 / sqrt(eSigma$values[1]))
  } else {
    inv_sqrt_vals <- 1 / sqrt(eSigma$values[1:r])
    halfCovInv <- eSigma$vectors[, 1:r, drop = FALSE] %*% diag(inv_sqrt_vals, nrow = r)
  }
  
  halfStat <- t(halfCovInv) %*% est
  stat <- sum(halfStat^2)
  p_val <- stats::pchisq(as.numeric(stat), df = r, lower.tail = FALSE)
  
  # 7. Calculate the AUC between fit_A and fit_B — LATE segment
  exp1 <- fit_A$prediction
  exp2 <- fit_B$prediction
  grid <- seq(0, 1, length.out = length(exp1))
  late_indices <- which(grid >= unname(interior_knots[n_knots]))  
  grid_late <- grid[late_indices]
  exp1_late <- exp1[late_indices]
  exp2_late <- exp2[late_indices]
  
  auc_late_A <- trapz(grid_late, exp1_late)
  auc_late_B <- trapz(grid_late, exp2_late)
  auc_difference <- log2(auc_late_A / auc_late_B)
  pred_val_A <- auc_late_A
  pred_val_B <- auc_late_B
  
  return(c(waldStat = as.numeric(stat), df = r, p_val = p_val, log2FC = auc_difference, predictA = pred_val_A, predictB = pred_val_B))
}

calculate_separate_wald <- function(fit_A, fit_B, l2fc= 0.4, eigenThresh = 1e-8) {
  # 1. Extract betas and vcovs
  betaA <- fit_A$beta
  betaB <- fit_B$beta
  VA <- fit_A$vcov
  VB <- fit_B$vcov
  L <- fit_A$X 
  
  # 2. Calculate the difference (Contrast)
  estFC <- L %*% betaA - L %*% betaB
  logFCCutoff <- log(2^l2fc) # log2 to log scale
  est <- sign(estFC)*pmax(0, abs(estFC) - logFCCutoff) # zero or remainder
  
  # 3. Calculate Combined Sigma (Independence assumption)
  # Sigma = L*VA*L' + L*VB*L'
  # Optimized as: L %*% (VA + VB) %*% t(L)
  sigma <- L %*% (VA + VB) %*% t(L)
  
  # 4. Eigen-decomposition
  eSigma <- eigen(sigma, symmetric = TRUE)
  
  # 5. Determine rank using the PatternTest threshold
  r <- try(sum(eSigma$values / eSigma$values[1] > eigenThresh), silent = TRUE)
  
  if (inherits(r, "try-error") || is.na(r) || r < 1) {
    return(c(waldStat = NA, df = NA, p_val = NA, log2FC = NA, predictA = NA, predictB = NA))
  }
  
  if (r == 1) {
    # If rank is 1, halfCovInv is just a column vector
    # sigma = v * lambda * v' => sigma^-1 = v * (1/lambda) * v'
    # halfCovInv = v * 1/sqrt(lambda)
    halfCovInv <- eSigma$vectors[, 1, drop = FALSE] * (1 / sqrt(eSigma$values[1]))
  } else {
    # Standard matrix logic for r > 1
    inv_sqrt_vals <- 1 / sqrt(eSigma$values[1:r])
    halfCovInv <- eSigma$vectors[, 1:r, drop = FALSE] %*% diag(inv_sqrt_vals, nrow = r)
  }
  # 6. Project the difference into the eigen-space
  # t(halfCovInv) is (r x N), estFC is (N x 1) -> halfStat is (r x 1)
  halfStat <- t(halfCovInv) %*% est
  stat <- sum(halfStat^2)
  p_val   <- stats::pchisq(as.numeric(stat), df = r, lower.tail = FALSE)
  # 7. Calculate the AUC between fit_A and fit_B
  exp1 = fit_A$prediction
  exp2 = fit_B$prediction
  grid = seq(0, 1, length.out = length(exp1))
  auc_dataset1 <- trapz(grid, exp1)
  auc_dataset2 <- trapz(grid, exp2)
  auc_difference <- log2(auc_dataset1/auc_dataset2)
  pred_val_A <- auc_dataset1
  pred_val_B <- auc_dataset2
  return(c(waldStat = as.numeric(stat), df = r, p_val = p_val, log2FC = auc_difference, predictA = pred_val_A, predictB = pred_val_B))
}

calculate_separate_diff <- function(fit_A, fit_B, l2fc = 0.1, eigenThresh = 1e-2) {
  # 1. Extract parameters
  betaA <- fit_A$beta
  betaB <- fit_B$beta
  VA <- fit_A$vcov
  VB <- fit_B$vcov
  
  # 2. L is the contrast
  L <- matrix(as.numeric(fit_A$X_end), ncol = 1)
  
  # 3. Combined Variance-Covariance Matrix
  Sigma <- VA + VB
  
  # 4. Project Difference and Variance into Contrast Space
  # estFC is the actual log-natural difference
  beta_diff <- betaA - betaB
  estFC <- t(L) %*% beta_diff
  
  # 5. sigma is the variance of that difference
  sigma <- t(L) %*% Sigma %*% L
  
  # 6. Apply the LFC Threshold (The tradeSeq "FC" logic)
  logFCCutoff <- log(2^l2fc)
  est <- sign(estFC) * pmax(0, abs(estFC) - logFCCutoff)
  
  # 7. Eigendecomposition Method (tradeSeq Stability Logic)
  eSigma <- eigen(sigma, symmetric = TRUE)
  
  # 8. Determine rank based on eigenvalue threshold
  # We use the ratio of eigenvalues to the maximum eigenvalue
  r <- sum(eSigma$values / eSigma$values[1] > eigenThresh)
  
  # 9. Handle cases where the matrix is numerically zero or invalid
  if (r == 0 || is.na(r)) {
    return(c(waldStat = NA, df = NA, p_val = NA, log2FC = NA, predictA = NA, predictB = NA))
  }
  
  # 10. Compute the pseudo-inverse component (halfCovInv)
  # This follows: halfCovInv = Vectors * diag(1/sqrt(Values))
  halfCovInv <- eSigma$vectors[, seq_len(r), drop = FALSE] %*% 
    diag(1 / sqrt(eSigma$values[seq_len(r)]), nrow = r)
  
  # 11. Final Wald Statistic
  # stat = t(est) %*% Sigma_Inv %*% est
  halfStat <- t(est) %*% halfCovInv
  wald_stat <- as.numeric(crossprod(t(halfStat)))
  
  # 12. Degrees of Freedom and P-value
  df_wald <- r
  p_val <- stats::pchisq(wald_stat, df = df_wald, lower.tail = FALSE)
  
  # 13. Convert natural log FC to Log2FC for reporting
  #log2FC_val <- as.numeric(estFC / log(2))
  
  pred_val_A <- fit_A$prediction[length(fit_A$prediction)]
  pred_val_B <- fit_B$prediction[length(fit_B$prediction)]
  log2FC_val <- log2(pred_val_A/pred_val_B)
  
  return(c(waldStat = wald_stat, df = df_wald, p_val = p_val, log2FC = log2FC_val, predictA = pred_val_A, predictB = pred_val_B))
}

run_pairwise_comparison <- function(list_A, list_B, interior_knots, n_knots) {
  # 1. Identify common genes
  common_genes <- intersect(names(list_A), names(list_B))
  if (length(common_genes) == 0) stop("No common genes found.")
  
  list_A <- list_A[common_genes]
  list_B <- list_B[common_genes]
  message(paste("Processing", length(common_genes), "common genes..."))
  
  # 2. Process genes in parallel
  results_list <- pbapply::pblapply(common_genes, function(gene) {
    fit_A <- list_A[[gene]]
    fit_B <- list_B[[gene]]
    
    # Template for failure (ensuring column count consistency)
    empty_res_p <- c(waldStat = NA, df = NA, p_val = NA, log2FC = NA, predictA = NA, predictB = NA)
    empty_res_d <- c(waldStat = NA, df = NA, p_val = NA, log2FC = NA, predictA = NA, predictB = NA)
    empty_res_d_part <- c(waldStat = NA, df = NA, p_val = NA, log2FC = NA, predictA = NA, predictB = NA)
    
    # --- Convergence and Null Checks ---
    # Check if beta or vcov exists
    total_failure <- is.null(fit_A$beta)  || is.null(fit_B$beta) || 
      is.null(fit_A$vcov)  || is.null(fit_B$vcov) ||
      isFALSE(fit_A$converged) || isFALSE(fit_B$converged)
    
    if (total_failure) {
      return(list(pattern = empty_res_p, diffend = empty_res_d, diffend_part = empty_res_d_part))
    }
    
    # --- Statistical Tests ---
    res_1 <- tryCatch({ 
      calculate_separate_wald(fit_A, fit_B) 
    }, error = function(e) empty_res_p)
    
    res_2 <- tryCatch({ 
      calculate_separate_diff(fit_A, fit_B) 
    }, error = function(e) empty_res_d)
    
    res_3 <- tryCatch({ 
      calculate_separate_wald_end(fit_A, fit_B, interior_knots = interior_knots, n_knots = n_knots) 
    }, error = function(e) empty_res_d)
    
    return(list(pattern = res_1, diffend = res_2, diffend_part = res_3))
  })
  
  # 3. Extract and format the Pattern Test Data Frame
  pattern_df <- as.data.frame(do.call(rbind, lapply(results_list, `[[`, "pattern")))
  rownames(pattern_df) <- common_genes
  colnames(pattern_df) <- c("waldStat", "df", "pvalue", "log2FC", "predictA", "predictB")
  
  # 4. Extract and format the DiffEnd Test Data Frame
  diffend_df <- as.data.frame(do.call(rbind, lapply(results_list, `[[`, "diffend")))
  rownames(diffend_df) <- common_genes
  colnames(diffend_df) <- c("waldStat", "df", "pvalue", "log2FC", "predictA", "predictB")
  
  # 5. Extract and format the DiffEndpart Test Data Frame
  diffend_part_df <- as.data.frame(do.call(rbind, lapply(results_list, `[[`, "diffend_part")))
  rownames(diffend_part_df) <- common_genes
  colnames(diffend_part_df) <- c("waldStat", "df", "pvalue", "log2FC", "predictA", "predictB")
  
  
  return(list(
    pattern = pattern_df,
    diffend = diffend_df,
    diffend_part = diffend_part_df
  ))
}

run_all_pairwise_bp <- function(cds, lineages, nknot = 8, n_knots = nknot-3) {
  
  pseudotime <- prepare_pt_flat(cds, lineages)
  interior_knots <- quantile(pseudotime, probs = seq(0, 1, length.out = nknot)[-c(1, nknot)])
  
  lineage_names <- lineages
  if (length(lineage_names) < 2) {
    stop("You need at least 2 lineages to perform pairwise comparisons.")
  }
  
  # 1. Generate all unique pairs of lineages
  pairs <- combn(lineage_names, 2, simplify = FALSE)
  
  # 2. Storage for results
  final_pattern_list <- list()
  final_diffend_list <- list()
  final_diffend_part_list <- list()
  
  for (pair in pairs) {
    l1 <- pair[1]
    l2 <- pair[2]
    
    # Get the sequence of Branch Point names for both lineages
    # (Assuming these are the names of the list elements in @expectation)
    branches_l1 <- names(cds@expectation[[l1]])
    branches_l2 <- names(cds@expectation[[l2]])
    
    # 3. Identify the shared path and the latest shared Branch Point (BP)
    shared_branches <- intersect(branches_l1, branches_l2)
    
    if (length(shared_branches) == 0) {
      warning(sprintf("No shared branches found between %s and %s. Skipping.", l1, l2))
      next
    }
    
    # The 'latest' shared branch is the last one in the common sequence
    latest_shared_bp <- tail(shared_branches, 1)
    comp_name <- paste(l1, "vs", l2, "at", latest_shared_bp, sep = "_")
    
    message(sprintf("\n>>> Comparing: %s vs %s", l1, l2))
    message(sprintf(">>> Common Branch Point identified: %s", latest_shared_bp))
    
    # 4. Extract data for that specific Branch Point from both lineages
    # We navigate into the expectation list using the lineage name, then the BP name
    data_l1 <- cds@expectation[[l1]][[latest_shared_bp]]
    data_l2 <- cds@expectation[[l2]][[latest_shared_bp]]
    
    # 5. Run existing comparison logic
    # Ensure run_pairwise_comparison handles the data objects found at these branch levels
    res <- run_pairwise_comparison(data_l1, data_l2, interior_knots = interior_knots, n_knots = n_knots)
    
    # 6. Store results using the specific comparison name
    final_pattern_list[[comp_name]] <- res$pattern
    final_diffend_list[[comp_name]] <- res$diffend
    final_diffend_part_list[[comp_name]] <- res$diffend_part
  }
  
  return(list(patterns = final_pattern_list, diffend = final_diffend_list, diffend_part = final_diffend_part_list))
}

reorder_metrics_grouped <- function(df, type) {
  cn <- colnames(df)
  
  # 1. Identify global columns
  global_cols <- c("waldStat", "pvalue", "df")
  global_present <- global_cols[global_cols %in% cn]
  
  # Identify pairwise columns (everything else)
  pairwise <- setdiff(cn, global_present)
  
  # 2. Extract base comparison names
  base_names <- sub("^(waldStat_|pvalue_|df_|log2FC_|predictA_|predictB_)", "", pairwise)
  unique_comps <- unique(base_names)
  
  # 3. Natural sort
  if (requireNamespace("gtools", quietly = TRUE)) {
    comp_order <- gtools::mixedsort(unique_comps)
  } else {
    comp_order <- sort(unique_comps)
  }
  
  # 4. Reconstruct columns in grouped order
  wald_cols      <- paste0("waldStat_", comp_order)
  pval_cols      <- paste0("pvalue_",   comp_order)
  df_cols        <- paste0("df_",       comp_order)
  l2fc_cols      <- paste0("log2FC_",   comp_order)
  predict_A_cols <- paste0("predictA_",   comp_order)
  predict_B_cols <- paste0("predictB_",   comp_order)
  
  # Filter for columns that actually exist
  ordered_pairwise <- c(wald_cols, pval_cols, df_cols, l2fc_cols, predict_A_cols, predict_B_cols)
  ordered_pairwise <- ordered_pairwise[ordered_pairwise %in% cn]
  
  # 5. Final Order
  final_order <- c(global_present, ordered_pairwise)
  df <- df[, final_order, drop = FALSE]
  cn_fin <- colnames(df)
  
  # 6. Dynamic Suffix Naming based on 'type'
  is_global <- cn_fin %in% c("waldStat", "pvalue", "df")
  
  # Dynamically adds '_combined_[type]' or '_[type]' to column names
  colnames(df) <- ifelse(is_global, 
                         paste0(cn_fin, "_combined_", type), 
                         paste0(cn_fin, "_", type))
  
  # Clean up formatting (_vs_ to vs)
  colnames(df) <- gsub("_vs_", "vs", colnames(df))
  
  return(df)
}


quasi_test_bp <- function(cds){
  lineages <- names(cds@lineages)
  
  # Both 'res' and 'waldResOmnibus' are now lists with names 'pattern' and 'diffend'
  res <- run_all_pairwise_bp(cds, lineages)
  
  for (lin in lineages) {
    # 1. Extract pairwise dataframes for this lineage
    # Assuming 'extract_lineage_df' handles one test type at a time
    out_pattern <- extract_lineage_df_bp(res$pattern, lin)
    out_diffend <- extract_lineage_df_bp(res$diffend, lin)
    out_diffend_part <- extract_lineage_df_bp(res$diffend_part, lin)
    
    cds@lineage_genes[[lin]] <- list() 
    cds@lineage_genes[[lin]]$lineage_genes <- list()
    
    # 1. Define the data frames and their corresponding types
    datasets <- list(
      pattern = out_pattern,
      diffend = out_diffend,
      diffend_part = out_diffend_part
    )
    
    # 2. Loop through and process each one
    processed_results <- list()
    
    for (type_name in names(datasets)) {
      # Run your custom reorder function
      df_test <- reorder_metrics_grouped(datasets[[type_name]], type = type_name)
      
      # Dynamically search for the log2FC columns for this specific type
      fc_regex <- sprintf("^log2FC_.*_%s$", type_name)
      fc_cols  <- grep(fc_regex, colnames(df_test), value = TRUE)
      
      # Dynamically define the new output column name (e.g., average_FC_start)
      target_col <- paste0("average_FC_", type_name)
      
      # Calculate the average log2FC
      if (length(fc_cols) == 1) {
        df_test[[target_col]] <- df_test[[fc_cols]]
      } else if (length(fc_cols) > 1) {
        df_test[[target_col]] <- apply(df_test[, fc_cols, drop = FALSE], 1, get_average_FC)
      } else {
        df_test[[target_col]] <- NA # Fallback if no log2FC columns are found
      }
      
      # Save back to the results list
      processed_results[[paste0(type_name, "_test")]] <- df_test
    }
    
    # 3. (Optional) Extract them back to individual global variables if needed
    pattern_test = processed_results$pattern_test
    diffend_test = processed_results$diffend_test
    diffend_part_test = processed_results$diffend_part_test
    
    
    cds@lineage_genes[[lin]]$lineage_genes$quasipoisson <- list(
      pattern_test = pattern_test,
      diffend_test = diffend_test,
      diffend_part_test = diffend_part_test
    )
  }
  return(cds)
}

extract_lineage_df_bp <- function(res, lineage) {
  nm <- names(res)
  
  # 1. Identify which list elements contain the lineage of interest
  keep <- sapply(nm, function(x) {
    # Split by "_at_" first to remove branch info, then by "_vs_"
    lineage_part <- strsplit(x, "_at_")[[1]][1]
    parts <- strsplit(lineage_part, "_vs_")[[1]]
    lineage %in% parts
  })
  
  sub_res <- res[keep]
  if (length(sub_res) == 0) return(NULL)
  
  out_list <- list()
  for (i in seq_along(sub_res)) {
    full_name <- names(sub_res)[i]
    df        <- sub_res[[i]]
    
    # 2. Robust Splitting: Separate Lineages and Branch
    main_parts   <- strsplit(full_name, "_at_")[[1]]
    lineage_part <- main_parts[1]
    branch_info  <- main_parts[2] # e.g., "aligned_BP_Y_189"
    
    lineage_pair <- strsplit(lineage_part, "_vs_")[[1]]
    a <- lineage_pair[1] # Reference
    b <- lineage_pair[2] # Comparison
    
    # 3. Standardization Logic
    # If the user-requested lineage is 'b', flip the data so it becomes 'a'
    if (lineage == b) {
      # Flip log2FC if it exists
      if ("log2FC" %in% colnames(df)) {
        df$log2FC <- -df$log2FC
      }
      
      # Swap PredictA and PredictB
      if ("predictA" %in% colnames(df) && "predictB" %in% colnames(df)) {
        temp_A <- df$predictA
        df$predictA <- df$predictB
        df$predictB <- temp_A
      }
      other <- a
    } else {
      other <- b
    }
    
    # 4. Create a descriptive new name that includes the branch
    # Format: Lineage_vs_Other_at_Branch
    new_name <- paste0(lineage, "_vs_", other, "_at_", branch_info)
    
    # Rename columns to include the specific comparison and branch info
    colnames(df) <- paste0(colnames(df), "_", new_name)
    out_list[[new_name]] <- df
  }
  
  # 5. Combine results
  # We use the full names to ensure columns are unique across different branches
  names(out_list) <- NULL
  out <- do.call(cbind, out_list)
  return(out)
}

get_average_FC <- function(FC){
  linear_values <- 2^FC
  mean_linear <- mean(linear_values, na.rm = TRUE)
  mean_log2 <- log2(mean_linear)
  mean_log2
}

get_meta_p <- function(Ps) {
  Ps <- Ps[!is.na(Ps)]
  Ps[Ps == 0] <- .Machine$double.xmin
  if (length(Ps) == 0) return(NA_real_)
  
  stat <- -2 * sum(log(Ps))
  df <- 2 * length(Ps)
  pchisq(stat, df=df, lower.tail=FALSE)
}

prepare_pt_flat <- function(cds, lineages) {
  # This list will temporarily hold the vectors for each BP
  all_pts_collected <- list()
  
  for (lineage in lineages) {
    # Access the 'streched' list which contains the 5 BPs
    bp_list <- cds@pseudotime[[lineage]]$streched
    
    for (bp_name in names(bp_list)) {
      # 1. Get raw pseudotime (no scaling)
      pt_raw <- bp_list[[bp_name]]
      
      # 2. Create unique names to identify exactly where this value came from
      # Format: Lineage.BP.CellIndex (e.g., ExN1.BP1.42)
      names(pt_raw) <- paste(lineage, bp_name, seq_along(pt_raw), sep = ".")
      
      # 3. Store in the collector list using a unique key
      key <- paste(lineage, bp_name, sep = "_")
      all_pts_collected[[key]] <- pt_raw
    }
  }
  
  # Collapse all lineages and all BPs into one single long vector
  final_vector <- unlist(all_pts_collected, use.names = TRUE)
  
  return(final_vector)
}

fit.m3_3_wald <- function(exp.sel, pt, size_factor, predict_pt, model, N, interior_knots) {
  
  if (!requireNamespace("speedglm", quietly = TRUE)) {
    stop("speedglm not installed")
  }
  
  exp_data.sel <- data.frame(
    pseudotime  = as.numeric(pt),
    size_factor = as.numeric(size_factor),
    expression  = as.numeric(exp.sel)
  )
  
  # Try to fit the model and extract wald components
  result <- tryCatch({
    
    # 1. Fit the model
    fit_model <- speedglm::speedglm(
      model,
      data = exp_data.sel,
      family = quasipoisson(),
      acc = 1e-3,
      model = TRUE, # We need this temporarily to extract terms
      y = FALSE
    )
    
    # 2. Get Predicted Values (y-hat) for the 100 points
    grid_data <- data.frame(
      pseudotime  = as.numeric(predict_pt),
      size_factor = 1
    )
    pred <- predict(fit_model, newdata = grid_data, type = "response")
    
    # Determine number of points for Wald stats: 2 * nknots (number of coefficients)
    n_coefs <- length(coef(fit_model))
    # nPoints_wald <- 2 * (n_coefs - 1) # If nknots means just the spline bases
    nPoints_wald <- 2 * n_coefs # If nknots means total parameters as per the user logic
    
    # 2. Wald Grid (100 points for the statistics)
    # This is what you use for the p-value math
    grid_wald <- data.frame(
      pseudotime = seq(0, 1, length.out = nPoints_wald), 
      size_factor = 1
    )
    
    end_pt_value <- 1
    df_end <- data.frame(pseudotime = end_pt_value, size_factor = 1)
    
    end_pt_value_2 <- unname(interior_knots[5])
    df_end_2 <- data.frame(
      pseudotime = seq(end_pt_value_2, 1, length.out = nPoints_wald), 
      size_factor = 1
    )
    
    # 3. Extract Basis Matrix (X) for the 100 points
    # This ensures the Wald test uses the same basis as the fit
    model_terms <- terms(fit_model)
    clean_terms <- delete.response(model_terms)
    X_matrix <- model.matrix(clean_terms, data = grid_wald)
    X_end <- model.matrix(clean_terms, data = df_end)
    X_end_2 <- model.matrix(clean_terms, data = df_end_2)
    
    # 4. Return all pieces needed for Wald
    list(
      prediction = as.numeric(pred),
      beta       = coef(fit_model),
      vcov       = vcov(fit_model),
      X          = X_matrix,
      X_end      = X_end,
      X_end_2    = X_end_2,
      converged  = fit_model$convergence
    )
    
  }, error = function(e) {
    # Return NAs/NULLs if the gene fails to fit
    list(
      prediction = rep(NA_real_, N),
      beta       = NULL,
      vcov       = NULL,
      X          = NULL,
      X_end      = NULL,
      X_end_2    = NULL,
      converged  = FALSE
    )
  })
  
  return(result)
}

quasi_refit_by_branchpoint <- function(cds, cores = 1, N = 1000, nknot = 8) {
  lineages <- names(cds@lineages)
  pseudotime <- prepare_pt_flat(cds, lineages)
  interior_knots <- quantile(pseudotime, probs = seq(0, 1, length.out = nknot)[-c(1, nknot)])
  predict_pt <- seq(0, 1, length.out = N)
  
  for (lineage in lineages) {
    # Access the list of 5 branch point pseudotimes
    bp_list <- cds@pseudotime[[lineage]]$streched
    bp_names <- names(bp_list)
    
    message(paste("Processing lineage:", lineage, "with branch points:", paste(bp_names, collapse = ", ")))
    
    # Initialize the container for this lineage
    cds@expectation[[lineage]] <- list()
    #pseudotime <- prepare_pt_flat(cds, lineages)
    #interior_knots <- quantile(pseudotime, probs = seq(0, 1, length.out = nknot)[-c(1, nknot)])
    #predict_pt <- seq(0, 1, length.out = N)
    
    # Iterate through each Branch Point (BP) element
    for (bp in bp_names) {
      message(paste("  Fitting model for Branch Point:", bp))
      
      # 1. Extract expression and size factors
      meta_sum <- cds@expression[[lineage]][["sum"]]
      mat_m <- as.matrix(meta_sum[, 7:ncol(meta_sum)])
      size_factor <- meta_sum$size_factor
      genes <- colnames(mat_m)
      
      # 2. Get the specific pseudotime vector for this BP
      # Using the 'streched' values directly
      pt_vec <- bp_list[[bp]]
      
      # Handle potential NAs (cells not belonging to this specific BP path)
      valid_cells <- !is.na(pt_vec)
      d <- pt_vec[valid_cells]
      
      # Subset matrix and size factors to match valid cells
      mat_valid <- mat_m[valid_cells, ]
      size_factor <- size_factor[valid_cells]
      
      # 4. Define the formula
      # Note: 'pt' and 'size_factor' are passed as arguments to the fit function
      model <- as.formula(
        substitute(
          expression ~ splines::ns(pseudotime, knots = k, Boundary.knots = c(0, 1)) + offset(log(size_factor)),
          list(k = interior_knots)
        )
      )
      
      # 5. Parallel Fitting across genes
      fit_list <- pbapply::pblapply(setNames(seq_along(genes), genes), function(i) {
        fit.m3_3_wald(
          exp.sel = mat_valid[, i],
          pt = d,
          size_factor = size_factor,
          predict_pt = predict_pt,
          model = model,
          N = N,
          interior_knots = interior_knots
        )
      }, cl = cores)
      # 6. Save results into the branch point list
      cds@expectation[[lineage]][[bp]] <- fit_list
    }
  }
  return(cds)
}

within_lineage_DE_par <- function(cds, lineages = names(cds@lineages), parallel = F, BPPARAM = BPPARAM, test = "Moran", cores = 1){
  if(test == "tradeSeq"){
    out <- lapply(lineages, within_lineage_DE_trade, cds = cds, parallel = T, BPPARAM = BPPARAM)
  }
  else(
    out <- lapply(lineages, within_lineage_DE_Moran, cds = cds, cores = cores)
  )
  names(out) <- lineages
  cds@dynamic_genes <- out
  cds
}

make_time_nb <- function(n, k = 5) {
  stopifnot(n >= 2, k >= 1)
  k <- min(k, n - 1)
  nb <- vector("list", n)
  for (i in seq_len(n)) {
    lo <- max(1, i - k)
    hi <- min(n, i + k)
    nb[[i]] <- setdiff(lo:hi, i)
  }
  class(nb) <- "nb"
  attr(nb, "region.id") <- as.character(seq_len(n))
  nb
}

within_lineage_DE_Moran <- function(lineage, #name of the lineage to analyze
                                    cds, #metatracker object
                                    k = 5,
                                    cores = 1
){
  print(paste0("Testing lineage ", lineage))
  d =cds@expression[[lineage]][["mean"]]
  expr = t(as.matrix(sapply(d[,6:ncol(d)], as.numeric))) #a matrix of expression values, with genes in rows and cells in columns
  expr = expr[rownames(cds),]
  n_time <- ncol(expr)
  nb <- make_time_nb(n_time, k = k)
  lw <- nb2listw(nb, style = "W", zero.policy = TRUE)
  expr[expr < 0] <- 0
  keep_rows <- apply(expr, 1, function(x) {
    all(is.finite(x)) && var(x) > 0})
  expr_filtered <- expr[keep_rows, , drop = FALSE]
  res <- pblapply(seq_len(nrow(expr_filtered)), function(i) {
    mt <- moran.test(expr_filtered[i, ], lw, zero.policy = TRUE)
    c(I = unname(mt$estimate[["Moran I statistic"]]),
      p = mt$p.value)
  }, cl = cores)
  res <- as.data.frame(do.call(rbind, res))
  res$padj <- p.adjust(res$p, method = "fdr")
  rownames(res) <- rownames(expr_filtered)
  res <- res[order(res$padj, -res$I),]
  res
}

within_lineage_DE_trade <- function(lineage, #name of the lineage to analyze
                                    cds, #metatracker object
                                    conditions = NULL, #a vector of condition information
                                    nknots = 3, #number of knots used to fit the GAM
                                    pairwise=TRUE, #pairwise comparison between different conditions
                                    contrast_type = "end",
                                    parallel = F,
                                    BPPARAM = F
){
  print(paste0("Testing lineage ", lineage))
  d =cds@expression[[lineage]]
  #conditions = d[,conditions]
  counts = t(as.matrix(sapply(d[,4:ncol(d)], as.numeric))) #a matrix of expression values, with genes in rows and cells in columns
  counts = counts[rownames(cds),]
  pseudotime = as.matrix(cds@pseudotime[[lineage]]) #a matrix of pseudotime values, each row for a cell
  cell_wt<-as.matrix(rep(1,ncol(counts)),ncol=1)
  #gamList<-tradeSeq::fitGAM(counts=counts,conditions=conditions,nknots=nknots,pseudotime=pseudotime,cellWeights=cell_wt,sce=FALSE)
  #res<-tradeSeq::conditionTest(gamList,global=TRUE,pairwise=pairwise,lineages=FALSE)
  #res_full<-res[!is.na(res$pvalue),]
  #res_sig<-res_full[res_full$pvalue<p,]
  #p_value <- getSmootherPvalues(gamList)
  #p_value = p_value[p_value<p,]
  #statLineage <- getSmootherTestStats(gamList)
  #statLineage = statLineage[names(p_value),]
  #final_res = cbind(statLineage, p_value)
  #final_res= final_res[order(final_res[,1], decreasing = T),]
  gamList<-tradeSeq::fitGAM(counts=counts,conditions=conditions,nknots=nknots,pseudotime=pseudotime,cellWeights=cell_wt, parallel = parallel, BPPARAM =BPPARAM)
  res = associationTest(gamList, contrastType = contrast_type)
  exp_95 = matrix(,nrow = nrow(res),ncol = 0)
  for(lin in names(cds@lineages)){
    exp_lin = cds@expectation[[lin]]
    exp_lin = exp_lin[,rownames(res)]
    exp_95_lin = apply(exp_lin, 2, function(x) as.numeric(quantile(x, 0.95, na.rm= TRUE)))
    exp_95 = cbind(exp_95, exp_95_lin)
  }
  rownames(exp_95) <- rownames(res)
  exp_95_max = apply(exp_95, 1, max, na.rm = TRUE)
  expectation = cds@expectation[[lineage]]
  expectation = expectation[,rownames(cds)]
  FC_factor = apply(expectation, 2, function(x) as.numeric(quantile(x, 0.95, na.rm= TRUE)))/exp_95_max
  scaled_FC = res$meanLogFC*FC_factor
  res$scaled_FC <- scaled_FC                      
  res = res_sig[with(res, order(pvalue, -scaled_FC)), ]
  res
}

filter_by_expression_lineage <- function(cds, lineage, mode = "number", N = 100, ratio = 0.01){
  data = counts(cds)
  cells = as.character(cds@lineages[[lineage]][['name']])
  data.sub = data[,cells]
  if(mode == "ratio"){
    cutoff = ratio*ncol(data.sub)
  }
  else{
    cutoff = N
  }
  expressed_genes = rownames(data.sub)[rowSums(data.sub> 0) > cutoff]
  expressed_genes
}


process_lineage <- function(df, adj_matrix, fc_matrix, p_cutoff, fc_cutoff, type = "quasipoisson") {
  
  # 1. Row filtering logic: 
  # Checks that every non-NA p-value is below cutoff AND every non-NA FC is above cutoff
  p_pass  <- rowSums(adj_matrix <= p_cutoff, na.rm = TRUE) == rowSums(!is.na(adj_matrix))
  fc_pass <- rowSums(fc_matrix >= fc_cutoff, na.rm = TRUE) == rowSums(!is.na(fc_matrix))
  
  # Combine conditions
  df <- df[p_pass & fc_pass, ]
  return(df)
}

format_quasipoisson_cols <- function(df) {
  # We use any_of() and matches() to ensure the function is robust 
  # even if some columns aren't present in every dataset.
  df_formatted <- df %>%
    select(
      matches("^waldStat|^log2FC"),
      starts_with("predictA"),
      starts_with("predictB"),
      starts_with("pvalue"),
      -any_of("pvalue_meta"),           # Temporarily remove to re-position at the end
      matches("^average_FC|^median_FC"),
      any_of("I"),
      any_of("pvalue_meta")             # Add it back at the very end
    )
  
  return(df_formatted)
}

format_lineage_specific_genes_bp <- function(lineage, cds, p_cutoff = 0.05, FC_pattern_cutoff = 0.2, FC_diffend_gen = 0.009, FC_diffend_cutoff = 0.4, dynamic_I_cutoff = 0.1, dynamic_p_cutoff = 0.05, p_adjust = "BH", specificity = "high", dynamic_test = "Moran", type = "quasipoisson"){
  lineages = names(cds@lineages)
  pattern_genes = cds@lineage_genes[[lineage]][["lineage_genes"]][[type]][['pattern_test']]
  diffend_genes = cds@lineage_genes[[lineage]][["lineage_genes"]][[type]][['diffend_test']]
  diffend_part_genes = cds@lineage_genes[[lineage]][["lineage_genes"]][[type]][['diffend_part_test']]
  #First filter out genes that express in fewer than 100 cells
  expressed_genes <- filter_by_expression_lineage(cds=cds, lineage=lineage, mode = "number", N = 100, ratio = 0.01)
  common_genes <- intersect(expressed_genes, rownames(pattern_genes))
  pattern_genes <- pattern_genes[common_genes, , drop = FALSE]
  common_genes <- intersect(expressed_genes, rownames(diffend_genes))
  diffend_genes <- diffend_genes[common_genes, , drop = FALSE]
  common_genes <- intersect(expressed_genes, rownames(diffend_part_genes))
  diffend_part_genes <- diffend_part_genes[common_genes, , drop = FALSE]
  #Second filter out genes that has p value >= 0.05 and moran's I statistics above 0.1
  if(dynamic_I_cutoff != F){
    dynamic = cds@dynamic_genes[[lineage]]
    if(dynamic_test == "Moran"){
      dynamic_genes = rownames(dynamic[dynamic$I >= dynamic_I_cutoff & dynamic$padj <= dynamic_p_cutoff, ])
    }
    else{
      dynamic_genes = rownames(dynamic[dynamic$scaled_FC >= dynamic_I_cutoff, ])
    }
    common_genes <- intersect(dynamic_genes, rownames(pattern_genes))
    pattern_genes <- pattern_genes[common_genes, , drop = FALSE]
    dynamic_subset_I <- dynamic[common_genes, "I", drop = FALSE]
    pattern_genes <- cbind(pattern_genes, I = dynamic_subset_I)
    
    common_genes <- intersect(dynamic_genes, rownames(diffend_genes))
    diffend_genes <- diffend_genes[common_genes, , drop = FALSE]
    dynamic_subset_I <- dynamic[common_genes, "I", drop = FALSE]
    diffend_genes <- cbind(diffend_genes, I = dynamic_subset_I)
    
    common_genes <- intersect(dynamic_genes, rownames(diffend_part_genes))
    diffend_part_genes <- diffend_part_genes[common_genes, , drop = FALSE]
    dynamic_subset_I <- dynamic[common_genes, "I", drop = FALSE]
    diffend_part_genes <- cbind(diffend_part_genes, I = dynamic_subset_I)
  }
  #Third filter out genes with NA values in averageFC
  pattern_genes_adjusted <- pattern_genes[!is.na(pattern_genes$average_FC_pattern),]
  diffend_genes_adjusted <- diffend_genes[!is.na(diffend_genes$average_FC_diffend),]
  diffend_part_genes_adjusted <- diffend_part_genes[!is.na(diffend_part_genes$average_FC_diffend_part),]
  
  #For diffend test, filter the endpoint values
  predictA_cols <- grep("^predictA_", colnames(diffend_genes_adjusted), value = TRUE)
  diffend_genes_adjusted <- diffend_genes_adjusted[rowSums(diffend_genes_adjusted[, predictA_cols, drop = FALSE] >= FC_diffend_gen, na.rm = TRUE) > 0, ]
  
  #predictA_cols <- grep("^predictA_", colnames(pattern_part_genes_adjusted), value = TRUE)
  #pattern_part_genes_adjusted <- pattern_part_genes_adjusted[rowSums(pattern_part_genes_adjusted[, predictA_cols, drop = FALSE] >= 0.03, na.rm = TRUE) > 0, ]
  
  if(p_adjust != FALSE){
    pairwise_columns_p <- grep("^pvalue_.*vs.*_pattern$", colnames(pattern_genes_adjusted), value = TRUE)
    pattern_genes_adjusted$pvalue_meta  <- apply(pattern_genes_adjusted[, pairwise_columns_p, drop = FALSE], 1, get_meta_p)
    lineage_genes_p <- pattern_genes_adjusted
    #remove the NAs in any p_matrix_sig columns
    #lineage_genes_p <- lineage_genes_p[complete.cases(lineage_genes_p[, pairwise_columns_p]), ]
    lineage_genes_p[, pairwise_columns_p] <- lapply(lineage_genes_p[, pairwise_columns_p, drop = FALSE], p.adjust, method = "BH")
    lineage_genes_p$pvalue_meta <- p.adjust(lineage_genes_p$pvalue_meta, method = "BH")
    p_matrix_complete_adj <- lineage_genes_p[, pairwise_columns_p, drop = FALSE]
    
    pairwise_columns_d <- grep("^pvalue_.*vs.*_diffend$", colnames(diffend_genes_adjusted), value = TRUE)
    diffend_genes_adjusted$pvalue_meta  <- apply(diffend_genes_adjusted[, pairwise_columns_d, drop = FALSE], 1, get_meta_p)
    lineage_genes_d <- diffend_genes_adjusted
    #lineage_genes_d <- lineage_genes_d[complete.cases(lineage_genes_d[, pairwise_columns_d]), ]
    lineage_genes_d[, pairwise_columns_d] <- lapply(lineage_genes_d[, pairwise_columns_d, drop = FALSE], p.adjust, method = "BH")
    lineage_genes_d$pvalue_meta <- p.adjust(lineage_genes_d$pvalue_meta, method = "BH")
    d_matrix_adj <- lineage_genes_d[, pairwise_columns_d, drop = FALSE]
    
    pairwise_columns_d_part <- grep("^pvalue_.*vs.*_diffend_part$", colnames(diffend_part_genes_adjusted), value = TRUE)
    diffend_part_genes_adjusted$pvalue_meta  <- apply(diffend_part_genes_adjusted[, pairwise_columns_d_part, drop = FALSE], 1, get_meta_p)
    lineage_genes_d_part <- diffend_part_genes_adjusted
    #lineage_genes_d_part <- lineage_genes_d_part[complete.cases(lineage_genes_d_part[, pairwise_columns_d_part]), ]
    lineage_genes_d_part[, pairwise_columns_d_part] <- lapply(lineage_genes_d_part[, pairwise_columns_d_part, drop = FALSE], p.adjust, method = "BH")
    lineage_genes_d_part$pvalue_meta <- p.adjust(lineage_genes_d_part$pvalue_meta, method = "BH")
    d_part_matrix_adj <- lineage_genes_d_part[, pairwise_columns_d_part, drop = FALSE]
  }
  if(length(lineages) > 2){
    FCs_p = lineage_genes_p[, grepl("log2FC", colnames(lineage_genes_p))]
    median_FC_p <- apply(FCs_p, 1, median, na.rm = TRUE)
    median_FC_p = median_FC_p[rownames(lineage_genes_p)]
    lineage_genes_p$median_FC_pattern <- median_FC_p
    median_FC_p = median_FC_p[rownames(pattern_genes_adjusted)]
    pattern_genes_adjusted$median_FC_pattern <- median_FC_p
    
    FCs_d = lineage_genes_d[, grepl("log2FC", colnames(lineage_genes_d))]
    median_FC_d <- apply(FCs_d, 1, median, na.rm = TRUE)
    median_FC_d = median_FC_d[rownames(lineage_genes_d)]
    lineage_genes_d$median_FC_diffend <- median_FC_d
    median_FC_d = median_FC_d[rownames(diffend_genes_adjusted)]
    diffend_genes_adjusted$median_FC_diffend <- median_FC_d
    
    FCs_d_part = lineage_genes_d_part[, grepl("log2FC", colnames(lineage_genes_d_part))]
    median_FC_d_part = apply(FCs_d_part, 1, median, na.rm = TRUE)
    median_FC_d_part = median_FC_d_part[rownames(lineage_genes_d_part)]
    lineage_genes_d_part$median_FC_diffend_part <- median_FC_d_part
    median_FC_d_part = median_FC_d_part[rownames(diffend_part_genes_adjusted)]
    diffend_part_genes_adjusted$median_FC_diffend_part <- median_FC_d_part
    
    if(specificity == "high"){
      lineage_genes_p <- process_lineage(lineage_genes_p, p_matrix_complete_adj, FCs_p, p_cutoff = p_cutoff, fc_cutoff = 0.2, type = "quasipoisson")
      lineage_genes_d_part <- process_lineage(lineage_genes_d_part, d_part_matrix_adj, FCs_d_part, p_cutoff = p_cutoff, fc_cutoff = 0.2, type = "quasipoisson")
      lineage_genes_d <- process_lineage(lineage_genes_d, d_matrix_adj, FCs_d, p_cutoff = p_cutoff, fc_cutoff = 0.4, type = "quasipoisson")
    }
  }
  lineage_genes_p <- format_quasipoisson_cols(lineage_genes_p)
  lineage_genes_d_part <- format_quasipoisson_cols(lineage_genes_d_part)
  lineage_genes_d <- format_quasipoisson_cols(lineage_genes_d)
  gene = union(rownames(lineage_genes_p), union(rownames(lineage_genes_d), rownames(lineage_genes_d_part)))
  # 1. Merge the datasets
  df_list <- list(
    lineage_genes_p[, c("pvalue_meta", "I"), drop = FALSE], 
    lineage_genes_d[, c("pvalue_meta", "I"), drop = FALSE], 
    lineage_genes_d_part[, c("pvalue_meta", "I"), drop = FALSE]
  )
  
  merged_df <- Reduce(function(x, y) {
    # If 'Row.names' already exists in x, merge by that column name. Otherwise use 0.
    by_x <- if("Row.names" %in% colnames(x)) "Row.names" else 0
    merge(x, y, by.x = by_x, by.y = 0, all = TRUE)
  }, df_list)
  
  rownames(merged_df) <- merged_df$Row.names
  colnames(merged_df) <- c("Row.names", "pvalue_p", "I_p", "pvalue_d", "I_d", "pvalue_d_part", "I_d_part")
  merged_df$I <- pmax(merged_df$I_p, merged_df$I_d, merged_df$I_d_part, na.rm = TRUE)
  merged_df$Row.names <- NULL
  merged_df <- merged_df[, c("pvalue_p", "pvalue_d", "pvalue_d_part", "I")]
  merged_df <- merged_df[gene, ]
  
  p_df <- pattern_genes_adjusted[gene, c("median_FC_pattern"), drop = FALSE]
  rownames(p_df) <- gene
  d_df <- diffend_genes_adjusted[gene, c("median_FC_diffend"), drop = FALSE]
  rownames(d_df) <- gene
  d_part_df <- diffend_part_genes_adjusted[gene, c("median_FC_diffend_part"), drop = FALSE]
  rownames(d_part_df) <- gene
  p_df <- p_df[rownames(d_df), , drop = FALSE]
  d_part_df <- d_part_df[rownames(d_df), , drop = FALSE]
  merged <- cbind(p_df, d_df, d_part_df, merged_df)
  colnames(merged) <- c("median_FC_pattern", "median_FC_diffend", "median_FC_diffend_part", "pvalue_meta_patterntest", "pvalue_meta_diffendtest", "pvalue_meta_diffendparttest", "moranI")
  merged$transientScore <- 
    rank(-merged$median_FC_pattern, ties.method = "min")^2 + rank(-merged$median_FC_diffend, ties.method = "min")^2 + rank(-merged$median_FC_diffend_part, ties.method = "min")^2
  lineage_spec_genes <- merged[order(merged$transientScore), ]
  lineage_spec_genes$lineage <- rep(lineage, nrow(lineage_spec_genes))
  lineage_spec <- list("pattern_filtered" = lineage_genes_p, "pattern_prefiltered" = pattern_genes_adjusted, "diffend_filtered" = lineage_genes_d, "diffend_prefiltered" = diffend_genes_adjusted, "diffend_part_filtered" = lineage_genes_d_part, "diffend_part_prefiltered" = diffend_part_genes_adjusted, "combined" = lineage_spec_genes)
  lineage_spec
}

format_branch_specific_genes_bp <- function(branch_point, cds, branch_number = 1, bp_name, p_cutoff = 0.05, FC_cutoff = 0.2){
  lineages = names(cds@lineages)
  branches_1 = branch_point[[branch_number]]
  if(branch_number == 1){
    branches_2 = branch_point[[2]]
  }
  else{
    branches_2 = branch_point[[1]]
  }
  branch_genes = list()
  branch_gene_names = list()
  for(lineage in branches_1){
    lineage_genes = cds@lineage_genes[[lineage]][["filtered"]][["quasipoisson"]][["pattern_prefiltered"]]
    FC_names = c()
    p_names = c()
    predict_names = c()
    for(lin in branches_2){
      FC_names = c(FC_names, paste0("log2FC_", lineage, "vs", lin, "_"))
      p_names = c(p_names, paste0("pvalue_", lineage, "vs", lin, "_"))
      predict_names = c(predict_names, paste0("predictA_", lineage, "vs", lin, "_"))
    }
    search_pattern <- paste(c(FC_names, p_names, predict_names), collapse = "|")
    lineage_genes_filtered <- lineage_genes[, grep(search_pattern, colnames(lineage_genes))]
    has_fc <- rowSums(!is.na(lineage_genes_filtered[, grepl("log2FC", colnames(lineage_genes_filtered)), drop = FALSE])) > 0
    has_p  <- rowSums(!is.na(lineage_genes_filtered[, grepl("pvalue", colnames(lineage_genes_filtered)), drop = FALSE])) > 0
    lineage_genes_filtered <- lineage_genes_filtered[has_fc & has_p, ,drop = FALSE]
    fc_cols <- grep("log2FC", colnames(lineage_genes_filtered), value = TRUE)
    p_cols  <- grep("pvalue", colnames(lineage_genes_filtered), value = TRUE)
    predictA_cols  <- grep("predictA", colnames(lineage_genes_filtered), value = TRUE)
    
    #pass_predict <- rowMeans(lineage_genes_filtered[, predictA_cols, drop = FALSE] >= 0.01, na.rm = TRUE) == 1
    pass_predict <- rowSums(lineage_genes_filtered[, predictA_cols, drop = FALSE] >= 0.01, na.rm = TRUE) > 0
    lineage_genes_filtered <- lineage_genes_filtered[which(pass_predict), ,drop = FALSE]
    
    lineage_genes_filtered$pvalue_meta  <- apply(lineage_genes_filtered[, p_cols, drop = FALSE], 1, get_meta_p)
    lineage_genes_filtered[, p_cols] <- lapply(lineage_genes_filtered[, p_cols, drop = FALSE], p.adjust, method = "BH")
    lineage_genes_filtered$pvalue_meta <- p.adjust(lineage_genes_filtered$pvalue_meta, method = "BH")
    
    pass_fc <- rowMeans(lineage_genes_filtered[, fc_cols, drop = FALSE] >= FC_cutoff, na.rm = TRUE) == 1
    pass_p  <- rowMeans(lineage_genes_filtered[, p_cols, drop = FALSE] <= p_cutoff, na.rm = TRUE) == 1
    lineage_spec_genes <- lineage_genes_filtered[which(pass_fc & pass_p), ,drop = FALSE]
    
    branch_genes[[lineage]] <- lineage_spec_genes
    branch_gene_names[[lineage]] <- rownames(lineage_spec_genes)
  }
  branch_gene_names = Reduce(intersect, branch_gene_names)
  median_FC_matrix = matrix(0,nrow = length(branch_gene_names), ncol = 0)
  meta_p_matrix = matrix(0,nrow = length(branch_gene_names), ncol = 0)
  for(lineage in branches_1){
    lineage_genes = branch_genes[[lineage]]
    FC_names = c()
    
    for(lin in branches_2){
      FC_names = c(FC_names, paste0("log2FC_", lineage, "vs", lin, "_"))
    }
    search_fc_pattern <- paste(c(FC_names), collapse = "|")
    FCs <- lineage_genes[branch_gene_names, grep(search_fc_pattern, colnames(lineage_genes)), drop = FALSE]
    median_FC_matrix = cbind(median_FC_matrix, FCs)
    pvalue <- lineage_genes[branch_gene_names, 'pvalue_meta', drop = FALSE]
    meta_p_matrix = cbind(meta_p_matrix, pvalue)
    branch2_names_str <- paste(branches_2, collapse = "&")
    colnames(meta_p_matrix)[ncol(meta_p_matrix)] <- paste0("pvalue_meta_", lineage, " vs ", branch2_names_str)
  }
  if(ncol(median_FC_matrix)==1){
    median_FC <- setNames(median_FC_matrix[,1], rownames(median_FC_matrix))
  }else{
    median_FC = apply(median_FC_matrix, 1, median)
  }
  meta_p_matrix <- meta_p_matrix[names(median_FC), ,drop = FALSE]
  final_out = cbind(median_FC, meta_p_matrix)
  final_out = as.data.frame(final_out)
  rownames(final_out) <- names(median_FC)
  lineage_string <- paste(branches_1, collapse = "_")
  branch_label <- paste0(bp_name, "_", lineage_string)
  final_out$branch <- branch_label
  colnames(final_out) <- c("median_FC", colnames(meta_p_matrix), "specific_to_lineages")
  final_out = final_out[with(final_out, order(-abs(median_FC))), ]
  final_out
}
