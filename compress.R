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
