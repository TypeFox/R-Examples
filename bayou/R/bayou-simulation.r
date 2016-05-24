#' Simulates parameters from bayou models
#' 
#' \code{priorSim} Simulates parameters from the prior distribution specified by \code{make.prior}
#' 
#' @param prior A prior function created by \code{bayou::make.prior}
#' @param tree A tree of class 'phylo'
#' @param plot A logical indicating whether the simulated parameters should be plotted
#' @param nsim The number of parameter sets to be simulated
#' @param ... Parameters passed on to \code{plotSimmap(...)}
#' 
#' @return A list of bayou parameter lists
#' 
#' @export
priorSim <- function(prior,tree,plot=TRUE,nsim=1, ...){
  tree <- reorder(tree,'postorder')
  model <- attributes(prior)$model
  dists <- attributes(prior)$dist
  fixed <- which(attributes(prior)$dist=="fixed")
  notfixed <- which(attributes(prior)$dist!="fixed")
  dists <- dists[notfixed]
  prior.params <- attributes(prior)$param
  rdists <- lapply(dists,function(x) gsub('^[a-zA-Z]',"r",x))
  prior.params <- lapply(prior.params,function(x) x[-which(names(x)=="log")])
  rdists.fx <- lapply(rdists,get)
  rdists.fx <- lapply(1:length(rdists.fx),function(x) .set.defaults(rdists.fx[[x]],defaults=prior.params[[x]]))
  names(rdists.fx) <- gsub('^[a-zA-Z]',"r",names(rdists))
  N <- sapply(names(rdists.fx),function(x) switch(x, ralpha=nsim, rsig2=nsim, rsig2jump=nsim, rhalflife=nsim, rVy=nsim, rh2=nsim, rP=nsim, rw2=nsim, rNe=nsim, rk=nsim, rtheta=NULL, rloc=NULL, rsb=NULL))
  varN <- which(sapply(N,is.null))
  N <- N[-varN]
  simpar <- lapply(1:nsim,function(i){ y <- lapply(names(N), function(x) rdists.fx[[x]](1)); names(y) <- gsub('^[a-zA-Z]',"",names(N)); y})
  if("dalpha" %in% names(fixed)){
    simpar <- lapply(simpar, function(x){x$alpha <- 0; x})
  }
  if("dsig2" %in% names(fixed)){
    simpar <- lapply(simpar, function(x){x$sig2 <- 0; x})
  }
  #if(model %in% c("OUcpp","QGcpp","OUreparcpp")){
  #  T <- sum(tree$edge.length[!(1:length(tree$edge.length) %in% exclude.branches)])
  #  pp <- tree$edge.length/T
  #  pp[1:length(pp) %in% exclude.branches] <- 0
  #  k <- lapply(simpar,function(x) x$k)
  #  sb <- lapply(k,function(x) .sample(1:length(tree$edge.length),x,replace=TRUE,prob=pp))
  #  loc <- lapply(sb,function(x) runif(length(x),min=0,max=tree$edge.length[x]))
  #  t2 <- lapply(k,function(x) 2:(x+1))
  #  simpar <- lapply(1:nsim,function(x) c(simpar[[x]],list(sb=sb[[x]],loc=loc[[x]],t2=t2[[x]])))
  #  theta <- lapply(1:nsim,function(x) pars2simmap(simpar[[x]],tree,sim.theta=TRUE,root.theta=rdists.fx$rtheta(1))$pars$theta)
  #  simpar <- lapply(1:nsim, function(x) c(simpar[[x]],list(theta=theta[[x]],ntheta=length(theta[[x]]))))
  #}
  #if(model %in% c("OU","QG","OUrepar")){
  if("dk" %in% names(fixed)){
    k <- rep(0,length(simpar))
    simpar <- lapply(simpar, function(x){x$k <- 0; x})
    sb <- lapply(k, function(x) numeric(x))
  } else {
    k <- sapply(simpar,function(x) x$k)
    sb <- lapply(k,function(x) rdists.fx$rsb(x))
  }
  if("dloc" %in% names(fixed)){
    loc <- lapply(1:length(k), function(x) 0.5*tree$edge.length[sb[[x]]])
  } else {
    loc <- lapply(1:length(k),function(x) rdists.fx$rloc(k[x])*tree$edge.length[sb[[x]]])
  }
  if(k > 0){
    t2 <- lapply(k,function(x) 2:(x+1))
  } else {t2 <- lapply(1:length(k), function(x) numeric(0))}
  theta <- lapply(k,function(x) rdists.fx$rtheta(x+1))
  simpar <- lapply(1:nsim,function(x) c(simpar[[x]],list(ntheta=k[x]+1, theta=theta[[x]],sb=sb[[x]],loc=loc[[x]],t2=t2[[x]])))
  #}
  if(plot){
    if(nsim>1){
      par(ask=TRUE)
    }
    for(i in 1:nsim){
      maps <- pars2simmap(simpar[[i]],tree)
      col <- maps$col
      plotSimmap(maps$tree,colors=col, ...)
    }
    
  }
  return(list(pars=simpar,tree=tree))
}

#' Simulates data from bayou models
#' 
#' \code{dataSim} Simulates data for a given bayou model and parameter set
#' 
#' @param pars A bayou formated parameter list
#' @param model The type of model specified by the parameter list (either "OU", "OUrepar" or "QG").
#' @param tree A tree of class 'phylo'
#' @param map.type Either "pars" if the regimes are taken from the parameter list, or "simmap" if taken from the stored simmap in the tree
#' @param SE A single value or vector equal to the number of tips specifying the measurement error that should be simulated at the tips
#' @param phenogram A logical indicating whether or not the simulated data should be plotted as a phenogram
#' @param ... Optional parameters passed to \code{phenogram(...)}.
#'
#' @description This function simulates data for a given set of parameter values.
#' 
#' @export
dataSim <- function(pars, model, tree, map.type="pars", SE=0, 
                    phenogram=TRUE, ...){
  if(model %in% c("QG")){
    pars$alpha <- QG.alpha(pars)
    pars$sig2 <- QG.sig2(pars)
  }
  if(model %in% c("OUrepar")){
    p <- OU.repar(pars)
    pars$alpha <- p$alpha
    pars$sig2 <- p$sig2
  }
  if(map.type=="simmap"){
    print("Using mapped regimes from ape tree file")
    maps <- tree$maps
  }
  if(map.type=="pars"){
    print("Using mapped regimes from parameter list")
    maps <- pars2simmap(pars, tree)$tree$maps
  }
  dummy <- rep(0, length(tree$tip.label))
  names(dummy) <- tree$tip.label
  tree$maps <- maps
  cache <- .prepare.ou.univariate(tree,dummy)
  cache$maps <- maps
  W <- .simmap.W(cache,pars)
  if(pars$k > 0){
    E.th <- W%*%pars$theta
  } else {E.th <- W*pars$theta}
  Sigma <- .ouMatrix(vcv.phylo(tree),pars$alpha)*pars$sig2
  diag(Sigma) <- diag(Sigma)+SE
  X <- mvrnorm(1,E.th,Sigma)
  if(phenogram){
    col <- c(1,rainbow(pars$k))
    names(col) <- 1:pars$ntheta
    phenogram(cache$phy,X,colors=col,spread.labels=FALSE, ...)
  }
  return(list(W=W, E.th=E.th,dat=X))
}

