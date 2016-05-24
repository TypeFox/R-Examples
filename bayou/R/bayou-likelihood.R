#' bayOU internal function. 
#' 
#' \code{.fix.root.bm} is an internal function and not generally called by the user
#' 
#' This is an internal function from geiger.
.fix.root.bm <- function (root, cache) {
  rtidx = cache$root
  cache$y["m", rtidx] = root
  attr(cache$y, "given")[rtidx] = as.integer(TRUE)
  cache
}
#geiger:::.fix.root.bm

#' bayOU internal function. 
#' 
#' \code{.ou.cache} is an internal function and not generally called by the user
#' 
#' This is an internal function that modifies the internal function \code{.ou.cache} in geiger for efficiency.
.ou.cache.fast <- function (cache) 
{
  ht = cache$ht
  N = cache$n.tip
  Tmax = ht$start[N + 1]
  mm = match(1:nrow(ht), cache$edge[, 2])
  ht$t1 = Tmax - ht$end[cache$edge[mm, 1]]
  ht$t2 = ht$start - ht$end + ht$t1
  z = function(alpha) {
    if (alpha < 0) 
      stop("'alpha' must be positive valued")
    if (alpha == 0){
      bl = ht$t2-ht$t1
    } else {
      bl = (1/(2 * alpha)) * exp(-2 * alpha * (Tmax - ht$t2)) * 
        -(expm1(-2 * alpha * ht$t2)) - (1/(2 * alpha)) * 
        exp(-2 * alpha * (Tmax - ht$t1)) * -(expm1(-2 * 
                                                      alpha * ht$t1))
    }
    cache$len = bl
    cache
  }
  attr(z, "argn") = "alpha"
  return(z)
}

#' bayOU internal function. 
#' 
#' \code{fastbm.lik} is an internal function and not generally called by the user
#' 
#' This is an internal function that modifies the internal function \code{bm.lik} in geiger for efficiency.
.fastbm.lik <- function (cache, dat,SE = NA, model = "OU", ...) {
  cache$dat <- dat
  cache$y[1,][1:cache$ntips] <- dat
  #cache = .prepare.bm.univariate(phy, dat, SE = SE, ...)
  cache$ordering = attributes(cache$phy)$order
  cache$N = cache$n.tip
  cache$n = cache$n.node
  cache$nn = (cache$root + 1):(cache$N + cache$n)
  cache$intorder = as.integer(cache$order[-length(cache$order)])
  cache$tiporder = as.integer(1:cache$N)
  cache$z = length(cache$len)
  FUN = switch(model, OU = .ou.cache.fast(cache))
  ll.bm.direct = function(cache, sigsq, q = NULL, drift = NULL, 
                          se = NULL) {
    n.cache = cache
    given = attr(n.cache$y, "given")
    if (is.null(q)) {
      llf = FUN()
    }
    else {
      llf = FUN(q)
    }
    ll = llf$len
    dd = 0
    if (!is.null(drift)) 
      dd = drift
    adjvar = as.integer(attr(n.cache$y, "adjse"))
    adjSE = any(adjvar == 1)
    .xxSE = function(cache) {
      vv = cache$y["s", ]^2
      ff = function(x) {
        if (any(ww <- adjvar == 1)) {
          vv[which(ww)] = x^2
          return(vv)
        }
        else {
          return(vv)
        }
      }
      return(ff)
    }
    modSE = .xxSE(n.cache)
    vv = as.numeric(modSE(se))
    datC = list(len = as.numeric(ll), intorder = as.integer(n.cache$intorder), 
                tiporder = as.integer(n.cache$tiporder), root = as.integer(n.cache$root), 
                y = as.numeric(n.cache$y["m", ]), var = as.numeric(vv), 
                n = as.integer(n.cache$z), given = as.integer(given), 
                descRight = as.integer(n.cache$children[, 1]), descLeft = as.integer(n.cache$children[, 
                                                                                                      2]), drift = as.numeric(dd))
    parsC = as.numeric(rep(sigsq, n.cache$z))
    out = .Call("bm_direct", dat = datC, pars = parsC, PACKAGE = "bayou")
    loglik <- sum(out$lq)
    if (is.na(loglik)) 
      loglik = -Inf
    attr(loglik, "ROOT.MAX") = out$initM[datC$root]
    class(loglik) = c("glnL", class(loglik))
    return(loglik)
  }
  class(ll.bm.direct) <- c("bm", "dtlik", "function")
  fx_exporter = function() {
    attb = c()
    if (!is.null(qq <- argn(FUN))) {
      adjQ = TRUE
      attb = c(attb, qq)
    }
    else {
      adjQ = FALSE
    }
    attb = c(attb, "sigsq")
    if (any(attr(cache$y, "adjse") == 1)) {
      attb = c(attb, "SE")
    }
    if (model == "drift") {
      attb = c(attb, "drift")
    }
    cache$attb = attb
    lik <- function(pars, ...) {
      recache = function(nodes = NULL, root = 6, #ROOT.MAX=6
                         cache) {
        r.cache = cache
        if (root == 6) { #ROOT.MAX=6
          rtmx = TRUE
        }
        else if (root %in% c(3,4)) {#ROOT.OBS=3; ROOT.GIVEN=4
          rtmx = FALSE
          r.cache$attb = c(cache$attb, "z0")
        }
        else {
          stop("unusable 'root' type specified")
        }
        r.cache$ROOT.MAX = rtmx
        if (!is.null(nodes)) {
          m = r.cache$y["m", ]
          s = r.cache$y["s", ]
          g = attr(r.cache$y, "given")
          nn = r.cache$nn
          r.cache$y = .cache.y.nodes(m, s, g, nn, r.cache$phy, 
                                     nodes = nodes)
        }
        r.cache
      }
      rcache = recache(..., cache = cache)
      attb = rcache$attb
      if (missing(pars)) 
        stop(paste("The following 'pars' are expected:\n\t", 
                   paste(attb, collapse = "\n\t", sep = ""), sep = ""))
      pars = .repars(pars, attb)
      names(pars) = attb
      if (adjQ) 
        q = pars[[qq]]
      else q = NULL
      sigsq = pars[["sigsq"]]
      if ("SE" %in% attb) 
        se = pars[["SE"]]
      else se = NULL
      if ("drift" %in% attb) 
        drift = -pars[["drift"]]
      else drift = 0
      if ("z0" %in% attb) 
        rcache = .fix.root.bm(pars[["z0"]], rcache)
      ll = ll.bm.direct(cache = rcache, sigsq = sigsq, 
                        q = q, drift = drift, se = se)
      return(ll)
    }
    attr(lik, "argn") = attb
    attr(lik, "cache") <- cache
    class(lik) = c("bm", "function")
    lik
  }
  likfx = fx_exporter()
  return(likfx)
}

#' Function for calculating likelihood of an OU model in bayou using pruning algorithm or matrix inversion
#' 
#' @param pars A list of parameters to calculate the likelihood
#' @param tree A phylogenetic tree of class 'phylo'
#' @param X A named vector giving the tip data
#' @param SE A named vector or single number giving the standard errors of the data
#' @param model Parameterization of the OU model. Either "OU", "QG" or "OUrepar".
#' @param invert A logical indicating whether the likelihood should be solved by matrix inversion, rather than
#' the pruning algorithm. This is primarily present to test that calculation of the likelihood is correct. 
#' 
#' @details This function can be used for calculating single likelihoods using previously implemented methods. It is
#' likely to become deprecated and replaced by \code{bayou.lik} in the future, which is based on \code{phylolm}'s threepoint
#' algorithm, which works on non-ultrametric trees and is substantially faster.
#' 
#' @return A list returning the log likelihood ("loglik"), the weight matrix ("W"), the optima ("theta"),
#' the residuals ("resid") and the expected values ("Exp").
#' 
#' @export
OU.lik <- function(pars,tree,X,SE=0,model="OU", invert=FALSE){
  if(class(tree)=="phylo"){
    cache <- .prepare.ou.univariate(tree, X, SE=SE)
  } else {cache <- tree}
  dat <- cache$dat
  if(model=="QG"){
    pars$alpha <- QG.alpha(pars)
    pars$sig2 <- QG.sig2(pars)
  }
  if(model=="OUrepar"){
    repar <- OU.repar(pars)
    pars$alpha <- repar$alpha
    pars$sig2 <- repar$sig2
  }
  W <- .parmap.W(cache,pars)
  if(pars$ntheta>1){
    E.th <- W%*%pars$theta
  } else {E.th <- W*pars$theta}
  X.c<-dat-as.vector(E.th)
  if(invert){
    n <- cache$n
    V <- vcv(cache$phy)*pars$sig2
    ouV <- .ouMatrix(V, pars$alpha)
    diag(ouV) <- diag(ouV)+cache$SE^2
    detV <- determinant(ouV, logarithm=TRUE)
    loglik <- -0.5*(n*log(2*pi)+detV$modulus+t(X.c)%*%solve(ouV)%*%(X.c))[1,1]
    return(list(loglik=loglik,W=W,theta=pars$theta,resid=X.c,Exp=E.th))
  } else {
    lnL.fx<-.fastbm.lik(cache,X.c,SE=TRUE,model="OU")
  #lnL.fx<-bm.lik(cache$phy,X.c,SE=NA,model="OU")
    loglik <- lnL.fx(pars=c(pars$alpha,pars$sig2,0),root=4)#ROOT.GIVEN=4
    return(list(loglik=loglik,W=W,theta=pars$theta,resid=X.c,Exp=E.th))
  }
}

.OU.lik <- function(pars,cache,X,SE=0,model="OU"){
  if(model=="QG"){
    pars$alpha <- QG.alpha(pars)
    pars$sig2 <- QG.sig2(pars)
  }
  if(model=="OUrepar"){
    repar <- OU.repar(pars)
    pars$alpha <- repar$alpha
    pars$sig2 <- repar$sig2
  }
  W <- .parmap.W(cache,pars)
  if(pars$ntheta>1){
    E.th=W%*%pars$theta
  } else {E.th=W*pars$theta}
  X.c<-X-as.vector(E.th)
  lnL.fx<-.fastbm.lik(cache,X.c,SE=TRUE,model="OU")
  loglik <- lnL.fx(pars=c(pars$alpha,pars$sig2,0),root=4)#ROOT.GIVEN=4
  list(loglik=loglik,W=W,theta=pars$theta,resid=X.c,Exp=E.th)
}


#' Function for calculating likelihood of an OU model in bayou using the threepoint algorithm
#' 
#' @param pars A list of parameters to calculate the likelihood
#' @param cache A bayou cache object generated using .prepare.ou.univariate
#' @param X A named vector giving the tip data
#' @param model Parameterization of the OU model. Either "OU", "QG" or "OUrepar".
#' 
#' @details This function implements the algorithm of Ho and Ane (2014) implemented
#'  in the package \code{phylolm} for the \code{OUfixedRoot} model. It is faster 
#'  than the equivalent pruning algorithm in geiger, and can be used on non-
#'  ultrametric trees (unlike OU.lik, which is based on the pruning algorithm in
#'  geiger). 
#' @export
bayou.lik <- function(pars, cache, X, model="OU"){
  if(model=="QG"){
    pars$alpha <- QG.alpha(pars)
    pars$sig2 <- QG.sig2(pars)
  }
  if(model=="OUrepar"){
    repar <- OU.repar(pars)
    pars$alpha <- repar$alpha
    pars$sig2 <- repar$sig2
  }
  n <- cache$n
  X <- cache$dat
  #W <- C_weightmatrix(cache,pars)$W
  #if(pars$ntheta>1){
  #  E.th=W%*%pars$theta
  #} else {E.th=W*pars$theta}
  #X.c<-X-as.vector(E.th)
  X.c <- C_weightmatrix(cache, pars)$resid
  #tree.trans = .transf.branch.lengths(cache, model, pars$alpha)
  transf.phy <- C_transf_branch_lengths(cache, 1, X.c, pars$alpha)
  transf.phy$edge.length[cache$externalEdge] <- transf.phy$edge[cache$externalEdge] + cache$SE[cache$phy$edge[cache$externalEdge, 2]]^2*(2*pars$alpha)/pars$sig2
  #transf.phy$diagMatrix <- transf.phy$diagMatrix
  #P <- cbind(y)
  #comp4 <- .three.point.compute2(transf.phy, cache, P, Q=NULL)
  #comp3 <- three.point.compute(transf.phy$tree,P,Q=NULL,transf.phy$diagMatrix)
  #comp2 <- .three.point.compute(transf.phy, cache, P, Q=NULL)
  comp <- C_threepoint(list(n=n, N=cache$N, anc=cache$phy$edge[, 1], des=cache$phy$edge[, 2], diagMatrix=transf.phy$diagMatrix, P=X.c, root=transf.phy$root.edge, len=transf.phy$edge.length))
  #t(y) %*% inv.ouV %*% (y)*1/pars$sig2
  if(pars$alpha==0){
    inv.yVy <- comp$PP
    detV <- comp$logd
  } else {
    inv.yVy <- comp$PP*(2*pars$alpha)/(pars$sig2)
    detV <- comp$logd+n*log(pars$sig2/(2*pars$alpha))
  }
  #log(det(souV))  
  llh <- -0.5*(n*log(2*pi)+detV+inv.yVy)
  #return(list(loglik=llh, W=W,theta=pars$theta,resid=X.c,Exp=E.th, comp=comp, transf.phy=transf.phy))
  return(list(loglik=llh, theta=pars$theta,resid=X.c, comp=comp, transf.phy=transf.phy))
}

.pruningwise.branching.times <- function(phy, n, des, anc) {     
  xx <- numeric(phy$Nnode)
  interns <- which(phy$edge[, 2] > n)
  for (j in length(interns):1) {
    i = interns[j]
    xx[des[i] - n] <- xx[anc[i] - n] + phy$edge.length[i]
  }
  depth <- xx[phy$edge[1, 1] - n] + phy$edge.length[1]
  xx <- depth - xx
  names(xx) <- if (is.null(phy$node.label)) (n + 1):(n + phy$Nnode) else phy$node.label
  return(xx)
}

.pruningwise.distFromRoot <- function(phy, n, N) {   
  xx <- numeric(phy$Nnode+n)
  for (i in N:1) xx[phy$edge[i, 2]] <- xx[phy$edge[i, 1]] + phy$edge.length[i]
  names(xx) <- if (is.null(phy$node.label)) 1:(n + phy$Nnode) else phy$node.label
  return(xx)
}
