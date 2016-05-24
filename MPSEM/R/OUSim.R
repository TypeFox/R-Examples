#
### MPSEM package R functions and C wrappers.
#
EvolveOptimMarkovTree <- function(tp,tw,anc,p=1,root=tp$edge[1,1]) {
  nn <- length(tp$tip.label)+tp$Nnode
  if(nrow(tw) != ncol(tw))
    stop("Transition probability matrix (tw) must be a square matrix")
  if(anc > nrow(tw))
    stop("Ancestral state (anc) not defined in the transition probability matrix (tw).")
  if(any((abs(rowSums(tw)-1)) > sqrt(.Machine$double.eps)))
    warning("The sum of transition probabilities is not systematically 1.")
  if(root > nn)
    stop("Invalid parameter root.")
  res <- t(matrix(.C("EvolveQC",
                     as.integer(tp$edge[,1]),
                     as.integer(tp$edge[,2]),
                     as.integer(nrow(tp$edge)),
                     as.integer(nn),
                     nv = double(p*nn),
                     as.double(t(tw)),
                     as.integer(nrow(tw)),
                     as.integer(anc),
                     as.integer(p),
                     as.integer(root))$nv,p,nn))
  if(!is.null(tp$node.label)) {
    rownames(res) <- c(tp$tip.label,tp$node.label)
  } else {
    rownames(res) <- c(tp$tip.label,rep("",tp$Nnode))
  }
  colnames(res) <- paste("Trial",1:p,sep="_")
  return(res)
}
#
TraitOUsimTree <- function(tp,a,sigma,opt,p=1,root=tp$edge[1,1]) {
  nn <- length(tp$tip.label)+tp$Nnode
  if(root > nn)
    stop("Invalid parameter root.")
  if(length(opt) != nn)
    stop("Optima don't match the number of nodes.")
  res <- matrix(.C("OUsim",
                   as.integer(tp$edge[,1]),
                   as.integer(tp$edge[,2]),
                   as.integer(nrow(tp$edge)),
                   as.integer(nn),
                   as.double(tp$edge.length),
                   as.double(a[1]),
                   as.double(sigma[1]),
                   as.double(opt),
                   as.integer(p),
                   as.integer(root),
                   out=double(p*nn))$out,nn,p)
  if(!is.null(tp$node.label)) {
    rownames(res) <- c(tp$tip.label,tp$node.label)
  } else {
    rownames(res) <- c(tp$tip.label,rep("",tp$Nnode))
  }
  colnames(res) <- paste("Trial",1:p,sep="_")
  return(res)
}
#
OUvar <- function(d,a=0,theta=1,sigma=1) {
  nd <- length(d)
  w <- numeric(nd)
  a <- rep(a, length.out = nd)
  theta <- rep(theta, length.out = nd)
  sigma <- rep(sigma, length.out = nd)
  nz <- a != 0
  w[nz] <- (theta[nz]*(1-exp(-a[nz]*d[nz])))**2 + (sigma[nz]**2)*(1-exp(-2*a[nz]*d[nz]))/(2*a[nz])
  w[!nz] <- (sigma[!nz]**2)*d[!nz]
  return(w)
}
#
PEMvar <- function(d,a=0,psi=1) {
  nd <- length(d)
  a <- rep(a, length.out = nd)
  psi <- rep(psi, length.out = nd)
  return(.C("PEMvar",
            as.double(d),
            as.integer(nd),
            as.double(a),
            as.double(psi),
            res=double(nd))$res)
}
#
## Does not deal with semidefinite cases (use eigen instead of chol).
TraitVarGraphSim <- function(x,variance,distance="distance",p=1,...) {
  if(attr(x,"class") != "graph")
    stop("Parameter 'x' must be of class 'graph'")
  if(is.null(x$edge[[distance]]))
    stop("There is no property '",distance,"' for the edges of the graph.")
  if(is.null(x$vertex$species)) {
    B <- MPSEM::PEMInfluence(x,mroot = FALSE)
  } else {
    B <- MPSEM::PEMInfluence(x,mroot = FALSE)[x$vertex$species,]
  }
  n <- nrow(B)
  if(!missing(variance)) {
    fargs <- as.list(match.call())
    wfform <- formals(variance)
    wfform[[1L]] <- x$edge[[distance]]
    for(i in names(wfform)[-1L]) {
      if (!is.null(fargs[[i]])) {
        wfform[[i]] <- fargs[[i]]
      } else if(!is.null(x$edge[[i]])) {
        wfform[[i]] <- x$edge[[i]]
      }
    }
    return(matrix(rnorm(p*n,0,1),p,n)%*%chol(B%*%diag(do.call(what=variance,args=as.list(wfform)))%*%t(B)))
  } else {
    return(matrix(rnorm(p*n,0,1),p,n)%*%chol(B%*%t(B)))
  }
}
#
