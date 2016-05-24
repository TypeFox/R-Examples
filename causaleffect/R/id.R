id <-
function(y, x, P, G, to) {
  G.obs <- observed.graph(G)
  v <- get.vertex.attribute(G, "name")
  v <- to[which(to %in% v)]
  anc <- ancestors(y, G.obs, to)

  # line 1
  if (length(x) == 0) {
    P$sumset <- union(setdiff(v, y), P$sumset)
    P$var <- v
    return(P)
  }

  # line 2
  if (length(setdiff(v, anc)) != 0) {
    anc.graph <- induced.subgraph(G, anc)
    var.nonign <- getNonIgnorableNodes(P)
    P$sumset <- intersect(setdiff(v, anc), var.nonign)
    return(id(y, intersect(x, anc), P, anc.graph, to))
  }

  # line 3
  G.x.overbar <- subgraph.edges(G, E(G)[!to(x)], delete.vertices = FALSE)
  w <- setdiff(setdiff(v, x), ancestors(y, observed.graph(G.x.overbar), to))
  if (length(w) != 0) return(id(y, union(x, w), P, G, to))


  # line 4
  G.remove.x <- induced.subgraph(G, v[!(v %in% x)])
  s <- c.components(G.remove.x, to)
  if (length(s) > 1) {
    productlist <- list()
    for (i in 1:length(s)) productlist[[i]] <- id(s[[i]], setdiff(v, s[[i]]), P, G, to)
    return(probability(sumset = setdiff(v, union(y, x)), recursive = TRUE, children = productlist))
  } else {
    s <- s[[1]]
    
    # line 5 
    cG <- c.components(G, to)
    if (identical(cG[[1]], v)) {
      v.string <- paste(v, sep = "", collapse = ",")
      s.string <- paste(s, sep = "", collapse = ",")
      stop("Graph contains a hedge formed by C-forests of nodes: \n {", v.string , "} and {", s.string , "}.", call. = FALSE)
    }
   
    # line 6
    is.element <- FALSE
    for (i in 1:length(cG)) {
      if (identical(s, cG[[i]])) {
        is.element <- TRUE
        break
      }
    }
    if (is.element) {
      ind <- which(v %in% s)
      if (length(s) > 1) {
        productlist <- list()
        for (i in 1:length(s)) { 
          if (P$recursive) {
            productlist[[i]] <- parse.joint(P, s[i], v[0:(ind[i]-1)], v)
          } else {
            P.prod <- P
            P.prod$var <- s[i]
            P.prod$cond <- v[0:(ind[i]-1)]
            productlist[[i]] <- P.prod           
          }
        }  
        return(probability(sumset = setdiff(s, y), recursive = TRUE, children = productlist))
      } else {
        if (P$recursive) {
          P.prod <- parse.joint(P, s[1], v[0:(ind[1]-1)], v)
          P.prod$sumset <- union(P.prod$sumset, setdiff(s, y))
          return(P.prod)
        } else {
          P.prod <- P
          P.prod$var <- s[1]
          P.prod$cond <- v[0:(ind[1]-1)]
          P.prod$sumset <- union(P.prod$sumset, setdiff(s, y))
          return(P.prod)         
        }
      }
    }
    
    # line 7
    set.contain <- unlist(lapply(cG, FUN = function(x) setequal(intersect(x, s), s)))
    set.contain <- which(set.contain)[1]
    s <- cG[[set.contain]]
    productlist <- list()
    ind <- which(v %in% s)
    s.graph <- induced.subgraph(G, s)
    if (length(s) > 1) {
      for (i in 1:length(s)) { 
        P.prod <- P
        P.prod$var <- s[i]
        P.prod$cond <- v[0:(ind[i]-1)]
        productlist[[i]] <- P.prod
      }
      return(id(y, intersect(x, s), probability(recursive = TRUE, children = productlist), s.graph, to))
    } else {
      P.prod <- P
      P.prod$var <- s[1]
      P.prod$cond <- v[0:(ind[1]-1)]
      return(id(y, intersect(x, s), P.prod, s.graph, to))
    }
  }
}
