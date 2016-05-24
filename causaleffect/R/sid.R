sid <-
function(y, x, Pstar, D, to) {
  v.s <- get.vertex.attribute(D, "name")
  s <- v.s[which(vertex.attributes(D)$description == "S")]
  s <- to[which(to %in% s)]
  D.causal <- induced.subgraph(D, v.s[!(v.s %in% s)])
  v <- get.vertex.attribute(D.causal, "name")
  v <- to[which(to %in% v)]
  D.obs <- observed.graph(D.causal)
  D.s.obs <- observed.graph(D)
  anc <- ancestors(y, D.obs, to)
  anc.s <- ancestors(y, D.s.obs, to)

  # line 1
  if (length(x) == 0) {
    Pstar$sumset <- union(setdiff(v, y), Pstar$sumset)
    Pstar$var <- v
    return(Pstar)
  }

  # line 2
  if (length(setdiff(v, anc)) != 0) {
    anc.graph <- induced.subgraph(D, anc.s)
    var.nonign <- getNonIgnorableNodes(Pstar)
    Pstar$sumset <- intersect(setdiff(v, anc), var.nonign)
    return(sid(y, intersect(x, anc), Pstar, anc.graph, to))
  }

  # line 3
  D.x.overbar <- subgraph.edges(D.causal, E(D.causal)[!to(x)], delete.vertices = FALSE)
  w <- setdiff(setdiff(v, x), ancestors(y, observed.graph(D.x.overbar), to))
  if (length(w) != 0) return(sid(y, union(x, w), Pstar, D, to))

  # line 4
  D.remove.x <- induced.subgraph(D.causal, v[!(v %in% x)])
  cc <- c.components(D.remove.x, to)
  if (length(cc) > 1) {
    productlist <- list()
    for (i in 1:length(cc)) productlist[[i]] <- sid(cc[[i]], setdiff(v, cc[[i]]), Pstar, D, to)
    return(probability(sumset = setdiff(v, union(y, x)), recursive = TRUE, children = productlist))
  }

  # line 5  
  else {
    cc <- cc[[1]]
    cG.s <- sc.components(D, to)  
    cG <- lapply(cG.s, function(x) setdiff(x, s))

    # line 6
    A <- unobserved.graph(D)
    A <- subgraph.edges(A, E(A)[!to(x)], delete.vertices = FALSE)
    A <- as.matrix(get.adjacency(A))
    if (dSep(A, s, y, x)) return(probability(var = y, do = x))

    # line 7
    if (identical(cG[[1]], v)) {
      v.string <- paste(v, sep = "", collapse = ",")
      cc.string <- paste(cc, sep = "", collapse = ",")
      stop("Graph contains an s-hedge formed by s*-trees of nodes: \n {", v.string , "} and {", cc.string , "}.", call. = FALSE) 
    }

    # line 8
    is.element <- FALSE
    for (i in 1:length(cG)) {
      if (identical(cc, cG[[i]])) {
        is.element <- TRUE
        break
      }
    }
    if (is.element) {
      ind <- which(v %in% cc)
      if (length(cc) > 1) {
        productlist <- list()
        for (i in 1:length(cc)) { 
          if (Pstar$recursive) {
            productlist[[i]] <- parse.joint(Pstar, cc[i], v[0:(ind[i]-1)], v)
          } else {
            Pstar.prod <- Pstar
            Pstar.prod$var <- cc[i]
            Pstar.prod$cond <- v[0:(ind[i]-1)]
            productlist[[i]] <- Pstar.prod           
          }
        }  
        return(probability(sumset = setdiff(cc, y), recursive = TRUE, children = productlist))
      } else {
        if (Pstar$recursive) {
          Pstar.prod <- parse.joint(Pstar, cc[1], v[0:(ind[1]-1)], v)
          Pstar.prod$sumset <- union(Pstar.prod$sumset, setdiff(cc, y))
          return(Pstar.prod)
      } else {
          Pstar.prod <- Pstar
          Pstar.prod$var <- cc[1]
          Pstar.prod$cond <- v[0:(ind[1]-1)]
          Pstar.prod$sumset <- union(Pstar.prod$sumset, setdiff(cc, y))
          return(Pstar.prod)         
        }
      }
    }
    
    # line 9
    set.contain <- unlist(lapply(cG, FUN = function(x) setequal(intersect(x, cc), cc)))
    set.contain <- which(set.contain)[1]
    cc <- cG[[set.contain]]
    cc.s <- cG.s[[set.contain]]
    productlist <- list()
    ind <- which(v %in% cc)
    cc.graph <- induced.subgraph(D, cc.s)
    if (length(cc) > 1) {
      for (i in 1:length(cc)) { 
        Pstar.prod <- Pstar
        Pstar.prod$var <- cc[i]
        Pstar.prod$cond <- v[0:(ind[i]-1)]
        productlist[[i]] <- Pstar.prod
      }
      return(sid(y, intersect(x, cc), probability(recursive = TRUE, children = productlist), cc.graph, to))
    } else {
      Pstar.prod <- Pstar
      Pstar.prod$var <- cc[1]
      Pstar.prod$cond <- v[0:(ind[1]-1)]
      return(sid(y, intersect(x, cc), Pstar.prod, cc.graph, to))
    }
  }
}
