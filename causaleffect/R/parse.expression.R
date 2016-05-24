parse.expression <- 
function(P, to, G.Adj) {
  if (P$recursive) {
    parse.children <- sapply(P$children, FUN = function(x) (x$recursive | length(x$sumset) > 0))
    if (sum(parse.children) > 0) {
      for (i in which(parse.children)) {
        P$children[[i]] <- parse.expression(P$children[[i]], to, G.Adj)
      }
    }
    if (length(P$children) > 0) {
      parse.children <- sapply(P$children, FUN = function(x) (x$recursive | length(x$sumset) > 0))
       if (sum(parse.children) > 0) return(P)
    } else return(NULL)
  }

  if (length(P$sumset) == 0) return(P)

  if (!P$recursive) { 
    if (P$sumset == P$var) return(NULL)
    else return(P)
  }

  ord.children <- order(unlist(lapply(P$children, FUN = function(x) which(to == x$var))), decreasing = TRUE)
  ord.sum <- order(sapply(P$sumset, FUN = function(x) which(to == x)), decreasing = TRUE)
 
  P.sum <- P
  P.sum$children <- P.sum$children[ord.children]
  P.sum$sumset <- P.sum$sumset[ord.sum]
  P.sum <- simplify(P.sum, G.Adj)

  P.parse <- probability(recursive = TRUE, children = list())
  remove <- c()
  if (length(P.sum$sumset) > 0) {
    j <- 1
    for (i in 1:length(P.sum$children)) {
      if (length(intersect(P.sum$children[[i]]$var, P.sum$sumset)) == 0 & length(intersect(P.sum$children[[i]]$cond, P.sum$sumset)) == 0) {
        P.parse$children[[j]] <- P.sum$children[[i]]
        remove <- c(remove, i)
        j <- j + 1
     }
    }
  }

  P.sum$children[remove] <- NULL
  if (length(P.sum$children) > 0) P.parse$children[[length(P.parse$children) + 1]] <- P.sum
  if (length(P.parse$children) == 0) return(P.sum)
  return(P.parse)  
  return(P.sum)     
}

