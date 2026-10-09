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




