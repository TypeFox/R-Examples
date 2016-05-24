.verify_arguments <- function(mutmodel, locus, alleles) {

}

mutmodel_not_mut <- function(mutmodel, locus, alleles) {
  .verify_arguments(mutmodel, locus, alleles)

  dw <- mutmodel_dw_mut(mutmodel, locus, alleles)
  up <- mutmodel_up_mut(mutmodel, locus, alleles)
  
  not_mut <- 1 - dw - up

  return(not_mut)
}

mutmodel_dw_mut <- function(mutmodel, locus, alleles) {
  .verify_arguments(mutmodel, locus, alleles)
  
  pars <- mutmodel$mutpars[, locus]
  
  if (mutmodel$modeltype == 1L) {
    # SMM
    dw <- as.numeric(rep(pars[1L], length(alleles)))
    names(dw) <- alleles
    return(dw)
  } else if (mutmodel$modeltype == 2L) {
    # LMM
    dw <- as.numeric(unlist(lapply(alleles, function(a) pars[1L] / (1 + exp(pars[2L]*(pars[3L] - a))))))
    names(dw) <- alleles
    return(dw)
  } else if (mutmodel$modeltype == 3L) {
    # EMM
    dw <- as.numeric(unlist(lapply(alleles, function(a) 1/((1+exp(pars[1L]+pars[2L]*a))*(1+exp(pars[3L]+pars[4L]*a))))))
    names(dw) <- alleles
    return(dw)
  } else {
    stop("Unknown model type")
  }
  
  return(NULL)
}

mutmodel_up_mut <- function(mutmodel, locus, alleles) {
  .verify_arguments(mutmodel, locus, alleles)
  
  pars <- mutmodel$mutpars[, locus]
  
  if (mutmodel$modeltype == 1L) {
    # SMM
    up <- as.numeric(rep(pars[2L], length(alleles)))
    names(up) <- alleles
    return(up)
  } else if (mutmodel$modeltype == 2L) {
    # LMM    
    up <- as.numeric(unlist(lapply(alleles, function(a) pars[4L] / (1 + exp(pars[5L]*(pars[6L] - a))))))
    names(up) <- alleles
    return(up)
  } else if (mutmodel$modeltype == 3L) {
    # EMM
    up <- as.numeric(unlist(lapply(alleles, function(a) exp(pars[3L]+pars[4L]*a)/((1+exp(pars[1L]+pars[2L]*a))*(1+exp(pars[3L]+pars[4L]*a))))))
    names(up) <- alleles
    return(up)
  } else {
    stop("Unknown model type")
  }
  
  return(NULL)
}


approx_stationary_dist <- function(mutmodel, alleles) {
  .verify_arguments(mutmodel, 1L, alleles)
  
  Ps <- lapply(1L:ncol(mutmodel$mutpars), function(locus) {
    P <- matrix(0, ncol = length(alleles), nrow = length(alleles), dimnames = list(alleles, alleles))
    diag(P) <- mutmodel_not_mut(mutmodel, locus, alleles)
    for (i in 2L:ncol(P)) P[i-1L, i] <- mutmodel_dw_mut(mutmodel, locus, alleles[i])
    for (i in 1L:(nrow(P)-1L)) P[i+1L, i] <- mutmodel_up_mut(mutmodel, locus, alleles[i])
    # Hack to fixed finite range:
    P[nrow(P)-1L, ncol(P)] <- 1L - P[nrow(P), ncol(P)]
    P[2L, 1L] <- 1 - P[1L, 1L]
    
    return(P)
  })

  lapply(Ps, colSums)

  statdist <- function(P) {
    n <- nrow(P)
    p <- c(rep(0, n-1), 1)
    C <- matrix(0, nrow = n, ncol = n)
    diag(C) <- 1
    C <- C - P
    C[n, ] <- 1
    Q <- solve(C, p)
    
    names(Q) <- rownames(P)
    
    return(Q) 
  }

  statdists <- lapply(Ps, statdist)
  statdists <- do.call(rbind, statdists)
  statdists <- t(statdists)
  colnames(statdists) <- colnames(mutmodel$mutpars)
  
  return(statdists)
}


