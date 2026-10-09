plot_combined_graph <- function(cds, reduction_method = "UMAP", color_cells_by = NULL, alpha = 0.4,
                                label_clusters = FALSE, highlight_clusters = NULL,
                                other_color = "grey85", show_others = TRUE,
                                highlight_lineage = NULL, highlight_color = "orange",
                                highlight_size = 0.3, highlight_alpha = 0.6,
                                highlight_on_top = TRUE, graph_which = "subgraph") {
  
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
