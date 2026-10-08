# S4 class definition for the meta-tracker data set.

#' @export
setClass("metatracker_data_set", contains = "cell_data_set", slots=c(dynamic_genes="list", lineage_genes="list", graphs = "list", lineages="list", expression="list", expectation="list", pseudotime="list")) -> metatracker_data_set
