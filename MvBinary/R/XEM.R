OneXemBlock <- function(x, weight, alpha, tol){
  db <- length(alpha)
  n <- sum(weight)
  # Initialization
  delta <- sample(0:1, db, replace = TRUE)
  epsilon <- runif(db)
  # Computation of the posterior probabilities: beta corresponds to the bounds of the integral where the likelihood is constant
  beta <- alpha * delta + (1-alpha) * (1-delta)
  ord <- order(beta, decreasing = FALSE)
  beta <- beta[ord]
  # fb sont les valeurs des constantes pour chaque individus
  matSup <- matrix(0, db , db+1)
  matSup[upper.tri(matSup)] <- 1
  matInf <- 1-matSup
  # Proba per individuals
  partinf <- alpha[ord] * (1-epsilon[ord]) + epsilon[ord] * delta[ord] 
  partsup <- alpha[ord] * (1-epsilon[ord]) + epsilon[ord] * (1-delta[ord])   
  baseSup <- log(sweep(x[,ord], 2, partsup, "*") + sweep(1-x[,ord], 2, 1-partsup, "*"))
  baseInf <- log(sweep(x[,ord], 2, partinf, "*") + sweep(1-x[,ord], 2, 1 - partinf, "*"))
  fij <- exp(baseSup %*% matSup + baseInf%*% matInf)
  repere <- (c(beta, 1)-c(0, beta))
  probaInd <- fij %*% repere
  # Loglikelihood computation
  loglikeprec <- -Inf
  loglikeactu <- sum(log(probaInd)*weight)
  while ((loglikeactu - loglikeprec)>tol){
    # E step
    tik <- sweep( sweep(fij,2,repere,"*"), 1, probaInd, "/")  
    # M step
    for (j in 1:db){
      if (j==1) tij <- tik[,1] else tij <- rowSums(tik[,1:j])
      # notations:
      a <- sum(x[,ord[j]]*weight*tij) #                a sum_{i=1}^n w_ix_{ij}t_{ij} ou t_{ij} proba d'avoir l'indicatrice alpha_j<u_j
      b <- sum(x[,ord[j]]*weight*(1-tij)) #            b sum_{i=1}^n w_ix_{ij}(1-t_{ij})
      c <- sum((1-x[,ord[j]])*weight*tij)#             c sum_{i=1}^n w_i(1-x_{ij})(t_{ij})
      d <- sum((1-x[,ord[j]])*weight*(1-tij))#         d sum_{i=1}^n w_i(1-x_{ij})(1-t_{ij})
      al <- alpha[ord[j]]
      resdelta1 <- optimize(function(eps) a*log((1-eps)*al + eps) + b*log((1-eps)*al) + c*log((1-eps)*(1-al)) + d*log(1-al + eps*al), interval = c(0, 1), maximum=TRUE)
      resdelta0 <- optimize(function(eps) b*log((1-eps)*al + eps) + a*log((1-eps)*al) + d*log((1-eps)*(1-al)) + c*log(1-al + eps*al), interval = c(0, 1), maximum=TRUE)
      if (resdelta0$objective > resdelta1$objective){
        delta[ord[j]] <- 0
        epsilon[ord[j]] <- resdelta0$maximum
      } else {
        delta[ord[j]] <- 1 
        epsilon[ord[j]] <- resdelta1$maximum
      } 
    }    
    # Computation of the posterior probabilities:
    beta <- alpha * delta + (1-alpha) * (1-delta)
    ord <- order(beta, decreasing = FALSE)
    beta <- beta[ord]
    # Proba per individuals
    partinf <- alpha[ord] * (1-epsilon[ord]) + epsilon[ord] * delta[ord] 
    partsup <- alpha[ord] * (1-epsilon[ord]) + epsilon[ord] * (1-delta[ord])   
    baseSup <- log(sweep(x[,ord], 2, partsup, "*") + sweep(1-x[,ord], 2, 1-partsup, "*"))
    baseInf <- log(sweep(x[,ord], 2, partinf, "*") + sweep(1-x[,ord], 2, 1 - partinf, "*"))
    fij <- exp(baseSup %*% matSup + baseInf%*% matInf)
    repere <- (c(beta, 1)-c(0, beta))
    probaInd <- fij %*% repere
    # Loglikelihood computation
    loglikeprec <- loglikeactu
    loglikeactu <- sum(log(probaInd)*weight)
  }
#  if (loglikeprec>loglikeactu) warning("Error in the EM block algorithm")
  if (delta[1]==0) delta <-  1- delta
  return(list(epsilon=epsilon, delta=delta, loglike=loglikeactu))
}


XEMblock <- function(x, weight, alpha, tol, nbiter){
  allblock <- lapply(as.list(1:nbiter), function(it) OneXemBlock(x, weight, alpha, tol))
  loglike <- rep(NA, nbiter)
  for (it in 1:nbiter) loglike <- allblock[[it]]$loglike
  return(allblock[[which.max(loglike)]])
}


XEMmodel <- function(dataset, alpha, tol, nbinit.EM, model){
  delta <- epsilon <- rep(0, ncol(dataset))
  loglike <- 0
  nbparam <- length(alpha)
  for (b in unique(model)){
    if (sum(model==b)>1){
      tmp <- uniquecombs(dataset[,which(model==b)])
      weight <- as.numeric(table(attr(tmp,"index")))
      tmp <- XEMblock(as.matrix(tmp), weight, alpha[which(model==b)], tol, nbinit.EM)
      delta[which(model==b)] <- tmp$delta
      epsilon[which(model==b)] <- tmp$epsilon
      loglike <- loglike + tmp$loglike
      nbparam <- nbparam + sum(model==b)
    }else{
      delta[which(model==b)] <- 1
      epsilon[which(model==b)] <- 0
      loglike <- loglike + log(alpha[which(model==b)])*sum(dataset[,which(model==b)]) + log(1-alpha[which(model==b)])*sum(1-dataset[,which(model==b)])
    }
  }
  return(list(delta=delta, epsilon=epsilon, alpha=alpha, loglike=loglike, nbparam=nbparam, bic=loglike - nbparam*0.5*log(nrow(dataset)), model=model))
}