toCytoscape <-
function (igraphobj) {
  
  # Extract graph attributes
  if (requireNamespace("igraph",quietly=TRUE)) {graph_attr = igraph::graph.attributes(igraphobj)} else {stop("igraph package must be installed to use this function")}
  
  # Extract nodes
  node_count = length(igraph::V(igraphobj))
  if('name' %in% igraph::list.vertex.attributes(igraphobj)) {
    igraph::V(igraphobj)$id <- igraph::V(igraphobj)$name
  } else {
    igraph::V(igraphobj)$id <-as.character(c(1:node_count))
  }
  
  
  nodes <- igraph::V(igraphobj)
  nds = list()
  
  v_attr = igraph::vertex.attributes(igraphobj)
  v_names = igraph::list.vertex.attributes(igraphobj)
  
  for(i in 1:node_count) {
#     node_attr = list()
#     
#     for(j in 1:length(v_attr)) {
#       node_attr[j] = v_attr[[j]][i]  
#     }
#     names(node_attr) = v_names
    nds[[i]] = list(data = mapAttributes(v_names, v_attr, i))
  }
  
  edges <- igraph::get.edgelist(igraphobj)
  edge_count = igraph::ecount(igraphobj)
  e_attr <- igraph::edge.attributes(igraphobj)
  e_names = igraph::list.edge.attributes(igraphobj)
  
  attr_exists = FALSE
  e_names_len = 0
  if(identical(e_names, character(0)) == FALSE) {
    attr_exists = TRUE
    e_names_len = length(e_names)
  }
  e_names_len <- length(e_names)
  
  eds = list()
  for(i in 1:edge_count) {
    st = list(source=toString(edges[i,1]), target=toString(edges[i,2]))
    
    # Extract attributes
    if(attr_exists) {
      eds[[i]] = list(data=c(st, mapAttributes(e_names, e_attr, i)))
    } else {
      eds[[i]] = list(data=st)
    }
  }

  el = list(nodes=nds, edges=eds)
  
  x <- list(data = graph_attr, elements = el)
  if (requireNamespace("RJSONIO",quietly=TRUE)) {return (RJSONIO::toJSON(x))} else {stop("RJSONIO package must be installed to use this function")}
}
