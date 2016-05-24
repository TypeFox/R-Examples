trz <-
function(y, x, P, Z, J, D, to) {
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
    P$sumset <- union(setdiff(v, y), P$sumset)
    P$var <- v
    return(P)
  }

  # line 2
  if (length(setdiff(v, anc)) != 0) {
    anc.graph <- induced.subgraph(D, anc.s)
    var.nonign <- getNonIgnorableNodes(P)
    P$sumset <- intersect(setdiff(v, anc), var.nonign)
    return(trz(y, intersect(x, anc), P, Z, J, anc.graph, to))
  }

  # line 3
  D.x.overbar <- subgraph.edges(D.causal, E(D.causal)[!to(x)], delete.vertices = FALSE)
  w <- setdiff(setdiff(v, x), ancestors(y, observed.graph(D.x.overbar), to))
  if (length(w) != 0) return(trz(y, union(x, w), P, Z, J, D, to))

  # line 4
  D.remove.x <- induced.subgraph(D.causal, v[!(v %in% x)])
  cc <- c.components(D.remove.x, to)
  if (length(cc) > 1) {
    productlist <- list()
    for (i in 1:length(cc)) productlist[[i]] <- trz(cc[[i]], setdiff(v, cc[[i]]), P, Z, J, D, to)
    return(probability(sumset = setdiff(v, union(y, x)), recursive = TRUE, children = productlist))
  }

  # line 5  
  else {
    cc <- cc[[1]]
    cG.s <- sc.components(D, to)  
    cG <- lapply(cG.s, function(x) setdiff(x, s))

    # line 6
    if (!identical(cG[[1]], v)) {

      # line 7
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
            if (P$recursive) {
              productlist[[i]] <- parse.joint(P, cc[i], setdiff(v[0:(ind[i]-1)], J), v)
            } else {
              P.prod <- P
              P.prod$var <- cc[i]
              P.prod$cond <- setdiff(v[0:(ind[i]-1)], J)
              productlist[[i]] <- P.prod           
            }
          }  
          return(probability(sumset = setdiff(cc, y), recursive = TRUE, children = productlist))
        } else {
          if (P$recursive) {          
            P.prod <- parse.joint(P, cc[1], setdiff(v[0:(ind[1]-1)], J), v)
            P.prod$sumset <- union(P.prod$sumset, setdiff(cc, y))
            return(P.prod)
        } else {
            P.prod <- P
            P.prod$var <- cc[1]
            P.prod$cond <- setdiff(v[0:(ind[1]-1)], J)
            P.prod$sumset <- union(P.prod$sumset, setdiff(cc, y))
            return(P.prod)         
          }
        }
      }
    
      # line 8
      set.contain <- unlist(lapply(cG, FUN = function(x) setequal(intersect(x, cc), cc)))
      set.contain <- which(set.contain)[1]
      cc <- cG[[set.contain]]
      cc.s <- cG.s[[set.contain]]
      productlist <- list()
      ind <- which(v %in% cc)
      cc.graph <- induced.subgraph(D, cc.s)
      if (length(cc) > 1) {
        for (i in 1:length(cc)) { 
          P.prod <- P
          P.prod$var <- cc[i]
          P.prod$cond <- setdiff(v[0:(ind[i]-1)], J)
          productlist[[i]] <- P.prod
        }
        return(trz(y, intersect(x, cc), probability(recursive = TRUE, children = productlist), Z, J, cc.graph, to))
      } else {
        P.prod <- P
        P.prod$var <- cc[1]
        P.prod$cond <- setdiff(v[0:(ind[1]-1)], J)
        return(trz(y, intersect(x, cc), P.prod, Z, J, cc.graph, to))
      }
    }
 
    # line 9
    else {
    
      # line 10
      A <- unobserved.graph(D)
      A <- subgraph.edges(A, E(A)[!to(x)], delete.vertices = FALSE)
      A <- as.matrix(get.adjacency(A))
      if (dSep(A, s, y, x) & (length(intersect(Z, x)) != 0)) {
        xcapz <- intersect(x, Z)
        D.remove.xcapz <- induced.subgraph(D, v[!(v %in% xcapz)])
        P$star <- FALSE
        P$do <- xcapz
        P <- activate.interventions(P, FALSE, xcapz)
        return(trz(y, setdiff(x, Z), P, setdiff(Z, x), xcapz, D.remove.xcapz, to))
      }
 
      # line 11
      else {
        v.string <- paste(v, sep = "", collapse = ",")
        cc.string <- paste(cc, sep = "", collapse = ",")
        stop("Graph contains a zs-hedge formed by s*-trees of nodes: \n {", v.string , "} and {", cc.string , "}.", call. = FALSE) 
      }
    }
  }  
}
