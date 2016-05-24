# 2014-09-01 CJS Fixed bug when logitP.fixed is fixed in first position
# 2014-09-01 CJS Needed to deal with different behaviour between OPENBugs and JAGS when the logitP parameters may be fixed.
#    OpenBUGS does NOT include the fixed logits in the returned MCMC sampler; JAGS does.
#    Consequently, the expansion of the the logitP from the samplers for the fixed logits only has to be done in OpenBugs and not JAGS
#    There was also a subtle bug in dealing with the multinomial distribution where the length of p (that had to be padded to deal with an OpenBugs problem
#    had to have the indicies explicitly stated.


PredictivePosterior.TSPNDENP <- function (n1,
                                          m2.expanded,
                                          u2,
                                          logitP.fixed,
                                          p,
                                          U,
                                          Theta,
                                          Delta.max,
					  engine) {
#  Generate Predictive Posterior Plot (Bayesian p-value) given the data
#  for a TimeStratified Petersen with Diagonal Elements and error
#    n1, m2, u2  = vectors of input data
#    p, U, Theta  = matrix of values (rows=number of posterior samples, columns=strata)
#                  These are returned from the call to OpenBugs/ JAGS
#

#select.m2 <- !is.na(m2)
#select.u2 <- !is.na(u2)

  #browser()

  s <- length(n1)
  t <- length(u2)

  ## Interleave p and logitP.fixed, ignoring extra p's added at end
  ## CJS - 2014-09-01 In JAGS, this is not needed because it passes the logitP matrix alread padded
  #                   the correct dimension. Not sure why this now did this?
  p.bkp <- p

  # 2014-09-01. Fixed a problem when fixed p is in first position and 1:0 doesn't work properly
  if(any(!is.na(logitP.fixed[1:t])) & tolower(engine)=="openbugs"  ){  # the second condition is for JAGS
    for(j in which(!is.na(logitP.fixed[1:t]))){
      if(j==1){ p <- cbind(expit(logitP.fixed[1]), p)}  # code below fails when j==1
      if(j> 1){ p <- cbind(p[,1:(j-1)],  expit(logitP.fixed[j]), p[,-(1:(j-1))]) }
    }
  }
# browser()



## Transform saved iterations for theta from vectors to full movement matrices
Theta.bkp <- Theta

Theta <- lapply(1:nrow(Theta.bkp),function(k){
  M <- Theta.bkp[k,,]

  tmp <- matrix(0,nrow=s,ncol=t)

  for(i in 1:length(n1))
    tmp[i,i:min(t,i+Delta.max)] <-
      M[i,1:min(t-i+1,Delta.max+1)]

  tmp
})

## Simulate data for each iteration
#  browser()
simData <- lapply(1:nrow(p),function(k) simTSPNDE(n1,U[k,],p[k,],Theta[[k]]))

#browser()
## Compute discrepancy measures
discrep <- t(sapply(1:nrow(p),function(k){

  ## 1) Observed vs expected values for recaptures of marked fish
  ## a) Observed data
  temp1.o <- sqrt(m2.expanded[,1:t]) - sqrt(n1 * t(t(Theta[[k]]) * p[k,1:t]))
  d1.m2.o <- sum(temp1.o^2,na.rm=TRUE)

  ## b) Simulate data
  temp1.s <- sqrt(simData[[k]]$m2[,1:t]) - sqrt(n1 * t(t(Theta[[k]]) * p[k,1:t]))
  d1.m2.s <- sum(temp1.s^2,na.rm=TRUE)

  ## 2) Observed vs expected values for captures of unmarked fish
  ## a) Observed data
  temp2.o <- sqrt(u2) - sqrt(U[k,] * p[k,1:t])
  d1.u2.o <- sum(temp2.o^2,na.rm=TRUE)

  ## b) Simulate data
  temp2.s <- sqrt(simData[[k]]$u2) - sqrt(U[k,] * p[k,1:t])
  d1.u2.s <- sum(temp2.s^2,na.rm=TRUE)

  ## 3) Deviance (-2*log-likelihood)
  ## a) Observed data
  #browser()
  d2.m2.o <- -2 * sum(sapply(1:length(n1),function(i){
    cellProbs <- Theta[[k]][i,] * p[k,1:t]  # 2014-09-01 need to ignore extra p's at end which were needed for OPENbugs quirk
    cellProbs <- c(cellProbs,1-sum(cellProbs))

    dmultinom(m2.expanded[i,],n1[i],cellProbs,log=TRUE)
  }))

  d2.u2.o <- -2 * sum(dbinom(u2,U[k,],p[k,1:t],log=TRUE))

  d2.o <- d2.m2.o + d2.u2.o

  ## b) Simulated data
  d2.m2.s <- -2 * sum(sapply(1:length(n1),function(i){
    cellProbs <- Theta[[k]][i,] * p[k,1:t]  # 2014-09-01 ditto to previous fix
    cellProbs <- c(cellProbs,1-sum(cellProbs))


    dmultinom(simData[[k]]$m2[i,],n1[i],cellProbs,log=TRUE)
  }))

  d2.u2.s <- -2 * sum(dbinom(simData[[k]]$u2,U[k,],p[k,1:t],log=TRUE))

  d2.s <- d2.m2.s + d2.u2.s

  c(d1.m2.o, d1.m2.s, d2.m2.o, d2.m2.s,
    d1.u2.o, d1.u2.s, d2.u2.o, d2.u2.s,
    d1.m2.o+d1.u2.o, d1.m2.s+d1.u2.s, d2.o, d2.s)
}))

discrep
}

simTSPNDE <- function(n1,U,p,Theta){
  ## Simulate data from the TSPNDE model conditional on values of n and U.

  s <- length(n1)
  t <- length(U)

  ## 1) Simulate matrix of recoveries
  m2 <- t(sapply(1:length(n1),function(i){
    cellProbs <- Theta[i,] * p[1:t]
    cellProbs <- c(cellProbs,1-sum(cellProbs))
    if( any(cellProbs < 0)){browser()}
    rmultinom(1,n1[i],cellProbs)[1:t]
  }))

  ## 2) Add number of individuals not recovered to last column of m2
  m2 <- cbind(m2,n1-apply(m2,1,sum))

  ## 3) Simulate captures of unmarked fish
  u2 <- rbinom(t,U,p)

  return(list(m2=m2,u2=u2))
}
