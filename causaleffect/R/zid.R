zid <-
function(y, x, Z, K, J, P, G, to) {
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
    return(zid(y, intersect(x, anc), Z, K, J, P, anc.graph, to))
  }

  # line 3
  xkj <- intersect(v, union(x, union(K, J)))
  G.xkj.overbar <- subgraph.edges(G, E(G)[!to(xkj)], delete.vertices = FALSE)
  zw <- intersect(setdiff(setdiff(v, xkj), ancestors(y, observed.graph(G.xkj.overbar), to)), Z)
  w <- setdiff(setdiff(setdiff(v, xkj), ancestors(y, observed.graph(G.xkj.overbar), to)), Z)
  if (length(union(zw, w)) != 0) {
    K <- union(K, zw)
    P$do <- union(P$do, union(K, J))
    G.k.overbar <- subgraph.edges(G, E(G)[!to(K)], delete.vertices = FALSE)
    return(zid(y, union(x, w), setdiff(Z, zw), K, J, P, G.k.overbar, to))
  }

  # line 4
  G.remove.xkj <- induced.subgraph(G, v[!(v %in% xkj)])
  s <- c.components(G.remove.xkj, to)
  if (length(s) > 1) {
    productlist <- list()
    for (i in 1:length(s)) {
      J.prod <- union(J, intersect(Z, setdiff(v, s[[i]])))
      P.prod <- P
      P.prod$do <- union(P$do, union(K, J.prod))
      G.j.overbar <- subgraph.edges(G, E(G)[!to(J.prod)], delete.vertices = FALSE)
      productlist[[i]] <- zid(s[[i]], setdiff(setdiff(v, s[[i]]), Z), setdiff(Z, setdiff(v, s[[i]])), K, J.prod, P.prod, G.j.overbar, to)    
    }
    return(probability(sumset = setdiff(v, union(y, union(x, K))), recursive = TRUE, children = productlist))
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
            productlist[[i]] <- parse.joint(P, s[i], setdiff(v[0:(ind[i]-1)], union(K, J)), v)
          } else {
            P.prod <- P
            P.prod$do <- union(K, J)
            P.prod$var <- s[i]
            P.prod$cond <- setdiff(v[0:(ind[i]-1)], union(K, J))
            productlist[[i]] <- P.prod           
          }
        }  
        return(probability(sumset = setdiff(s, y), recursive = TRUE, children = productlist))
      } else {
        if (P$recursive) {
          P.prod <- parse.joint(P, s[1], setdiff(v[0:(ind[1]-1)], union(K, J)), v)
          P.prod$sumset <- union(P.prod$sumset, setdiff(s, y))
          return(P.prod)
        } else {
          P.prod <- P
          P.prod$do <- union(K, J)
          P.prod$var <- s[1]
          P.prod$cond <- setdiff(v[0:(ind[1]-1)], union(K, J))
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
        P.prod$cond <- setdiff(v[0:(ind[i]-1)], union(K, J))
        productlist[[i]] <- P.prod
      }
      return(zid(y, intersect(x, s), Z, K, J, probability(recursive = TRUE, children = productlist), s.graph, to))
    } else {
      P.prod <- P
      P.prod$var <- s[1]
      P.prod$cond <- setdiff(v[0:(ind[1]-1)], union(I, J))
      return(zid(y, intersect(x, s), Z, K, J, P.prod, s.graph, to))
    }
  }
}
