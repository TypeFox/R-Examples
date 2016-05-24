
name.proto <- function(., envir = parent.frame()) {
   stopifnot(is.environment(.) || 
      (is.character(.) && is.environment(get(., envir))))
   if (is.environment(.)) {
      if (exists("..Name", ., inherits = FALSE)) .$..Name
      else {
         L <- unlist(eapply(envir, identical, .))
         if (any(L)) names(L[L])[1]
         else gsub("^.* |>$", "", capture.output(print.default(.))[[1]])
      }
   } else {
      e <- get(., envir)
      if (exists("..Name", e, inherits = FALSE)) e$..Name
      else .
   }
}
      
graph.proto <- function(e,
   g = new("graphNEL", edgemode = "directed"), child.to.parent = TRUE) {
   if (missing(e)) e <- 
      if (exists(".that")) get(".that") else parent.frame()
   if ( suppressWarnings(! require(graph) || ! require(Rgraphviz)) )
      stop("Error: packages graph and Rgraphviz must be available for loading")
   # add node if its not already in g
   addNode. <- function(x, g) if (x %in% nodes(g)) g else addNode(x, g)
   # add edge between nodes adding nodes too if not already in g
   addEdge. <- function(x, y, g) {
      g <- addNode.(x,g); g <- addNode.(y, g)
      addEdge(x, y, g, 1)
   }
   nn <- unlist(eapply(e, is.proto))
   for(x in names(nn[nn]))
      g <- if (child.to.parent)
         addEdge.(name.proto(x,e), name.proto(get(x,e)$parent.env(), e), g)
      else
         addEdge.(name.proto(get(x,e)$parent.env(), e), name.proto(x, e), g)
   g
}

### test
# a <- proto()
# b <- a$proto()
# g <- graph.proto()
# plot(g)
# g <- graph.proto(child.to.parent = FALSE) # change arrow heads
# plot(g)
# g <- graph.proto(new("graphNEL")) # undirected
# plot(g)
# g <- graph.proto()
# attrs <- list(node = list(fillcolor = "lightgreen"), 
#               edge = list(color = "cyan"),
#		graph = list(rankdir = "BT"))
# plot(graph.proto(), attrs) # specify plot attributes

