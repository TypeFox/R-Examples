#' @title Plot a network with species and interactions highlighted
#' @description
#' Plot differences between two networks
#'
#' @param n1 a network
#' @param n2 a second network
#' @param na color of items unique to network 1
#' @param nb color of items unique to network 2
#' @param ns color of shared items
#' @param ... additional arguments to be passed to plot
#'
#' @return Nothing
#'
#' @export
network_betaplot <- function(n1, n2, na='#2ca02c', nb='#1f77b4', ns='grey', ...){
   M <- metaweb(list(n1,n2))
   s1 <- igraph::V(n1)$name
   s2 <- igraph::V(n2)$name
   # VERTICES
   vert_color <- rep(ns, length(igraph::V(M)))
   names(vert_color) <- igraph::V(M)$name
   for(v in igraph::V(M)$name)
   {
      if ((v %in% s1) && !(v %in% s2)) vert_color[v] = na
      if (!(v %in% s1) && (v %in% s2)) vert_color[v] = nb
   }
   # EDGES
   make_readable_edgelist = function(x) {
     X = igraph::get.edgelist(x)
     return(paste(X[,1], X[,2]))
   }
   edge_color <- rep(ns, length(igraph::E(M)))
   em = make_readable_edgelist(M)
   names(edge_color) <- em
   e1 = make_readable_edgelist(n1)
   e2 = make_readable_edgelist(n2)
   for(ee in em)
   {
      if ((ee %in% e1) && !(ee %in% e2)) edge_color[ee] = na
      if (!(ee %in% e1) && (ee %in% e2)) edge_color[ee] = nb
   }
   igraph::plot.igraph(M, vertex.color = vert_color, edge.color = edge_color, ...)
}
