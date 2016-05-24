usid <-
function(y, x, Pstar, D, to) {
  d <- length(D)
  v.s <- lapply(D, function(x) get.vertex.attribute(x, "name"))
  s <- lapply(1:d, function(x) v.s[[x]][which(vertex.attributes(D[[x]])$description == "S")])
  s <- lapply(1:d, function(x) to[[x]][which(to[[x]] %in% s[[x]])])
  D.causal <- induced.subgraph(D[[1]], v.s[[1]][!(v.s[[1]] %in% s[[1]])])
  v <- get.vertex.attribute(D.causal, "name")
  v <- to[[1]][which(to[[1]] %in% v)]
  D.obs <- observed.graph(D.causal)
  D.s.obs <- lapply(D, function(x) observed.graph(x))
  anc <- ancestors(y, D.obs, to[[1]])
  anc.s <- lapply(1:d, function(x) ancestors(y, D.s.obs[[x]], to[[x]]))

  # line 1
  if (length(x) == 0) {
    Pstar$sumset <- union(setdiff(v, y), Pstar$sumset)
    Pstar$var <- v
    return(Pstar)
  }

  # line 2
  if (length(setdiff(v, anc)) != 0) {
    anc.graph <- lapply(1:d, function(x) induced.subgraph(D[[x]], anc.s[[x]]))
    var.nonign <- getNonIgnorableNodes(Pstar)
    Pstar$sumset <- intersect(setdiff(v, anc), var.nonign)
    return(usid(y, intersect(x, anc), Pstar, anc.graph, to))
  }

  # line 3
  D.x.overbar <- subgraph.edges(D.causal, E(D.causal)[!to(x)], delete.vertices = FALSE)
  w <- setdiff(setdiff(v, x), ancestors(y, observed.graph(D.x.overbar), to[[1]]))
  if (length(w) != 0) return(usid(y, union(x, w), Pstar, D, to))

  # line 4
  D.remove.x <- induced.subgraph(D.causal, v[!(v %in% x)])
  cc <- c.components(D.remove.x, to[[1]])
  if (length(cc) > 1) {
    productlist <- list()
    for (i in 1:length(cc)) productlist[[i]] <- usid(cc[[i]], setdiff(v, cc[[i]]), Pstar, D, to)
    return(probability(sumset = setdiff(v, union(y, x)), recursive = TRUE, children = productlist))
  }

  # line 5  
  else {
    cc <- cc[[1]]
    cG.s <- lapply(1:d, function(x) sc.components(D[[x]], to[[x]]))  
    cG <- lapply(cG.s[[1]], function(x) setdiff(x, s[[1]]))

    # line 6
    for (i in 1:d) {
      A <- unobserved.graph(D[[i]])
      A <- subgraph.edges(A, E(A)[!to(x)], delete.vertices = FALSE)
      A <- as.matrix(get.adjacency(A))
      if (dSep(A, s[[i]], y, x)) return(probability(var = y, do = x, domain = letters[i]))
    }

    # line 7
    if (identical(cG[[1]], v)) {
      v.string <- paste(v, sep = "", collapse = ",")
      cc.string <- paste(cc, sep = "", collapse = ",")
      stop("Graph contains a hedge formed by sC-forests of nodes: \n {", v.string , "} and {", cc.string , "}.", call. = FALSE) 
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
    cc.s <- lapply(cG.s, function(x) x[[set.contain]])
    productlist <- list()
    ind <- which(v %in% cc)
    cc.graph <- lapply(1:d, function(x) induced.subgraph(D[[x]], cc.s[[x]]))
    if (length(cc) > 1) {
      for (i in 1:length(cc)) { 
        Pstar.prod <- Pstar
        Pstar.prod$var <- cc[i]
        Pstar.prod$cond <- v[0:(ind[i]-1)]
        productlist[[i]] <- Pstar.prod
      }
      return(usid(y, intersect(x, cc), probability(recursive = TRUE, children = productlist), cc.graph, to))
    } else {
      Pstar.prod <- Pstar
      Pstar.prod$var <- cc[1]
      Pstar.prod$cond <- v[0:(ind[1]-1)]
      return(sid(y, intersect(x, cc), Pstar.prod, cc.graph, to))
    }
  }
}
