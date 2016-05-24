# log likelihood
log.l <- function(thetas,
                  alpha,
                  hsnet,
                  ltheta,
                  together = together,
                  possible.stats,
                  GERGM_Object = GERGM_Object) {
  #turn dataframe into matrix
  hsnet <- as.matrix(hsnet)
  if (nrow(hsnet) == 1) {
    theta <- thetas[1:nrow(hsnet)]
    z <- hsnet * (theta - ltheta)
  }
  if (nrow(hsnet) > 1) {
    theta <- thetas[1:ncol(hsnet)]
    #print(str(hsnet))
    z <- hsnet %*% (theta - ltheta)
  }

  if(GERGM_Object@is_correlation_network){
    #this will calculate the h statistics on the original network as desired
    temp <- h.corr(possible.stats,
                   alpha,
                   theta = theta,
                   together = together,
                   GERGM_Object)[1, ]
  }else{
    dir <- TRUE
    if(!GERGM_Object@directed_network){
      dir <- FALSE
    }
    num.nodes <- nrow(GERGM_Object@bounded.network)
    triples <- t(combn(1:num.nodes, 3))
    temp <- h2(net = GERGM_Object@bounded.network,
       triples = triples,
       statistics = GERGM_Object@stats_to_use,
       alphas = GERGM_Object@weights,
       together = together,
       directed = dir)
  }
  ret <- rbind(theta) %*% temp - max(z) - log(sum(exp(z - max(z))))
  #print(ret)
  return(ret)
}

llg <- function(par,
                alpha,
                theta,
                z,
                together = together,
                possible.stats = possible.stats,
                GERGM_Object = GERGM_Object) {
  # log likelihood for unbounded network with g function
  statistics <- GERGM_Object@stats_to_use
  alphas <- GERGM_Object@weights
  net <- GERGM_Object@network
  beta <- par[1:(length(par) - 1)]
  sig <- 0.01 + exp(par[length(par)])
  BZ <- 0
  for (i in 1:(dim(z)[3])) {
    BZ <- BZ + beta[i] * z[, , i]
  }
  transformation_type <- GERGM_Object@transformation_type
  if(transformation_type == "logcauchy"){
    net2 <- pst(log(net), BZ, sig, 1)
    last_term <- sum(log(dst(log(net[upper.tri(net)]), BZ[upper.tri(net)], sig, 1))) +
      sum(log(dst(log(net[lower.tri(net)]), BZ[lower.tri(net)], sig, 1)))
  }
  if( transformation_type == "cauchy"){
    net2 <- pst(net, BZ, sig, 1)
    last_term <- sum(log(dst(net[upper.tri(net)], BZ[upper.tri(net)], sig, 1))) +
      sum(log(dst(net[lower.tri(net)], BZ[lower.tri(net)], sig, 1)))
  }
  if(transformation_type == "lognormal"){
    net2 <- pst(log(net), BZ, sig, Inf)
    last_term <- sum(log(dst(log(net[upper.tri(net)]), BZ[upper.tri(net)], sig, Inf))) +
      sum(log(dst(log(net[lower.tri(net)]), BZ[lower.tri(net)], sig, Inf)))
  }
  if( transformation_type == "gaussian"){
    net2 <- pst(net, BZ, sig, Inf)
    last_term <- sum(log(dst(net[upper.tri(net)], BZ[upper.tri(net)], sig, Inf))) +
      sum(log(dst(net[lower.tri(net)], BZ[lower.tri(net)], sig, Inf)))
  }
  num.nodes <- nrow(net2)
  triples <- t(combn(1:num.nodes, 3))

  log.li <- rbind(theta) %*%
    h2(net2,
       triples = triples,
       statistics = statistics,
       alphas = alphas,
       together = together) +
    last_term
  return(as.numeric(log.li))
}

# maximum pseudo-likelihood estimates
mple <- function(net, statistics, directed, verbose = TRUE) {
  xy <- net2xy(net, statistics, directed = directed)
  x <- xy$x
  y <- xy$y
  est <- coef(lm(y ~ x - 1))
  ests <- NULL
  if(verbose){
    ests <- optim(par = est, pl, y = y, x = x, method = "BFGS",
                  hessian = TRUE,control = list(fnscale = -1, trace = 6))
  }else{
    ests <- optim(par = est, pl, y = y, x = x, method = "BFGS",
                  hessian = TRUE,control = list(fnscale = -1, trace = 0))
  }
  return(ests)
}

# ------------------------------------------------------------
## Functions for correlation matrices
#1: jacobian of transformation of correlation matrices to the [0,1] space
jacobian <- function(partials){
  corrs.1 <- diag(partials[-nrow(partials), -1])
  d <- nrow(partials)
  prod.1 <- prod((1 - corrs.1^2)^(d-2))
  prod.2 <- 1
  for(k in 2 : (d - 2)){
    for(i in 1 : (d - k)){
      prod.2 = prod.2*(1-(partials[i,i+k])^2)^(d-1-k)
    }
  }
  result <- 2*((prod.1^(d-2))*prod.2)^(0.5)
  return(result)
}

#pseudo-likelihood for the correlation matrix
pl.corr <- function(theta, y, x, Jacobian){
  return(sum(log(dtexp(y, x %*% theta))) + log(Jacobian))
}

#MPLE for correlation matrices
mple.corr <- function(net,
                      bounded.net,
                      statistics,
                      directed = FALSE,
                      verbose = TRUE){
  xy.full <- net2xy(net, statistics, directed = directed)
  x <- xy.full$x #x's are the change statistics associated with the unbounded network
  xy.bounded <- net2xy(bounded.net, statistics, directed = directed)
  y <- xy.bounded$y #y's are the edge weights from the bounded [0,1] network
  J <- jacobian(bounded.net)
  est <- coef(lm(y ~ x - 1))
  ests <- NULL
  if(verbose){
    ests <- optim(par = est, pl.corr, y = y, x = x, Jacobian = J, method = "BFGS",
                  hessian = TRUE, control = list(fnscale = -1, trace = 6))
  }else{
    ests <- optim(par = est, pl.corr, y = y, x = x, Jacobian = J, method = "BFGS",
                  hessian = TRUE, control = list(fnscale = -1, trace = 0))
  }
  return(ests)
}

