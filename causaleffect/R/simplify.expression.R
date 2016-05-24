simplify.expression <-
function(P, G, to) {
  from <- NULL
  if (P$recursive) {
    for (i in 1:length(P$children)) {
      P$children[[i]] <- simplify.expression(P$children[[i]], G, to)
    }
  } else {
    if (!is.null(P$do)) {
      j <- 0
      while (j < length(P$do)) {
        j <- j + 1
        G.x <- subgraph.edges(G, E(G)[!to(P$do[-j])], delete.vertices = FALSE)
        zw <- setdiff(P$do[j], ancestors(P$cond, G.x, to)) 
        G.xzw <- subgraph.edges(G, E(G)[!to(P$do[-j]) & !to(zw)], delete.vertices = FALSE)
        G.xzw <- as.matrix(get.adjacency(G.xzw)) 

        # Application of rule 3 of do-calculus
        if (dSep(G.xzw, P$var, P$do[j], c(P$do[-j], P$cond))) {
          P$do <- P$do[-j]
          j <- 0
        } else {
          G.xz <- subgraph.edges(G, E(G)[!to(P$do[-j]) & !from(P$do[j])], delete.vertices = FALSE)
          G.xz <- as.matrix(get.adjacency(G.xz))  
      
          # Application of rule 2 of do-calculus        
          if (dSep(G.xz, P$var, P$do[j], c(P$do[-j], P$cond))) {
            P$cond <- union(P$cond, P$do[j])
            P$do <- P$do[-j]
            j <- 0
          }
        }
      }
    }
    j <- 0
    while (j < length(P$cond)) {
      j <- j + 1
      G.x <- subgraph.edges(G, E(G)[!to(P$do)], delete.vertices = FALSE)
      G.x <- as.matrix(get.adjacency(G.x))

      # Application of rule 1 of do-calculus  
      if (dSep(G.x, P$var, P$cond[j], union(P$cond[-j], P$do))) {
        P$cond <- P$cond[-j]
        j <- 0
      }
    } 
  }
  if (P$fraction) P$divisor <- simplify.expression(P$divisor, G, to)
  return(P)
}
