## ----compIssue1----------------------------------------------------------
 p <- 1*10^-17
 x <- log(p/(1-p))
 log(1 - exp(x)/(1 - exp(x)))

## ----compIssue2----------------------------------------------------------
p <- 1*10^-17
x <- log(p/(1-p))
-log1p(exp(x))

## ----tab1SuccDir---------------------------------------------------------
library(SALTSampler)
TargetDir <- function(ycand, ycurrent, a = NULL, dat = NULL) {
  ycand <- Logit(ycand)
  ycurrent <- Logit(ycurrent)
  out <- sum((a - 1)*(LogPq(ycand)$logp - LogPq(ycurrent)$logp))
  return(out)
}

## ----succDirProp---------------------------------------------------------
PropStepDir <- function(y, s, guard) {
  p <- length(y)
  ynew <- rep(NA, p)
  #Sample first y value 
  ynew[1] <- rbeta(n = 1, shape1 = y[1]*s, shape2 = sum(y[2:p]*s))
  #Sample next p - 2  y values
  for (i in 2:(p-1)) {
    ynew[i] <- rbeta(n = 1, shape1 = y[i]*s, 
                     shape2 = sum(y[(i + 1):p]*s))*(1 - sum(ynew[1:(i - 1)]))
  }
  #Sample last y value
  ynew[p] <- 1 - sum(ynew[1:(p - 1)])
  #Calculate detailed balance term
  if (any(ynew < guard)) {
    dbt <- NA #Outside guard, always reject
  } else {
    dbt <- lgamma(sum(ynew*s)) - sum(lgamma(ynew*s)) + sum((ynew*s - 1)*log(y))
           - lgamma(sum(y*s))
    + sum(lgamma(y*s)) - sum((y*s - 1)*log(ynew)) 
  }
  #Return new logit-scaled point and corresponding dbt
  attr(ynew, 'dbt') <- dbt
  return(ynew)
}

## ----DirSampler----------------------------------------------------------
RunMhDir <- function(center, B, concentration, s, type, dat = NULL, guard) {
  #Redefining parameters, setting initial time, and making empty vectors and matrices
  zz = proc.time();
  p <- length(center)
  Y <- array(0, c(B, p)) 
  moveCount <- 0
  center <- center/sum(center)
  a <- concentration*center 
  ycurrent <- center
  #Run sampler
  for (i in 1:B) { 
    #Propose step
    ycand <- PropStepDir(ycurrent, s, guard)
    #Decide to accept or reject
    if (any(ycand < guard)) {
      move <- FALSE #Outside guard, always reject
    } else {
      move <- (log(runif(1)) < attr(ycand, 'dbt') + TargetDir(ycand, ycurrent, a, dat))
    }
    if (!is.na(move) & move) {
      ycurrent <- ycand
      moveCount <- moveCount + 1
    }
    #Store y for this iteration
    Y[i, ] <- ycurrent
  }
  #Timing
  runTime <- proc.time() - zz
  #Return results
  return(list(Y = Y, runTime = runTime, moveCount = moveCount, p = p, 
              center = center, B = B, concentration = concentration, s = s, 
              type = type, dat = NULL, a = a, moveCount = moveCount))
}

## ----succDirRun----------------------------------------------------------
succDir <- RunMhDir(center = c(1/3, 1/3, 1/3), B = 5e3, concentration = 3, 
                    s = 1.2, type = 'dirichlet', dat = NULL, guard = 0.01)

## ----tab1----------------------------------------------------------------
tab1 <- matrix(nrow = 6, ncol = 3)
tab1[, 1] <- c("0.01-0.02", "0.02-0.03", "0.03-0.05", "0.05-0.1", ">0.1", "Overall")
colnames(tab1) <- c("Minimum Coordinate", "Number of Samples", "Acceptance Rate")
minVal <- apply(succDir$Y[1:succDir$B - 1, ], 1, min)
g1 <- which(minVal < .02)
tab1[1, 2] <- length(g1)
tab1[1, 3] <- round(sum(succDir$Y[g1, ] != succDir$Y[g1 + 1, ])/length(g1), 5)
g2 <- which(minVal < .03 & minVal >= 0.02)
tab1[2, 2] <- length(g2)
tab1[2, 3] <- round(sum(succDir$Y[g2, ] != succDir$Y[g2 + 1, ])/length(g2), 5)
g3 <- which(minVal < .05 & minVal >= 0.03)
tab1[3, 2] <- length(g3)
tab1[3, 3] <- round(sum(succDir$Y[g3, ] != succDir$Y[g3 + 1, ])/length(g3), 5)
g4 <- which(minVal < .1 & minVal >= .05)
tab1[4, 2] <- length(g4)
tab1[4, 3] <- round(sum(succDir$Y[g4, ] != succDir$Y[g4 + 1, ])/length(g4), 5)
g5 <- which( minVal >= 0.1)
tab1[5, 2] <- length(g5)
tab1[5, 3] <- round(sum(succDir$Y[g5, ] != succDir$Y[g5 + 1, ])/length(g5), 5)
tab1[6, 2] <- succDir$B - 1
tab1[6, 3] <- round(succDir$moveCount/(succDir$B - 1), 5)
tab1

## ----resProp-------------------------------------------------------------
PropStepRP <- function(y, s) {
  p <- length(y)
  ynew <- rep(NA, p)
  #Set one coordinate aside at random
  hold <- sample(c(1:p), 1)
  use <- 1:p
  use <- use[use != hold]
  #Sample other p - 1 coordinates
  for (i in sample(use)) {
    ynew[i] <- rnorm(1, y[i], s)
  }
  #Calculate coordinate p from other coordinates
  ynew[hold] <- 1 - sum(ynew[use])
  #Calculate detailed balance term
  if (sum(ynew) != 1 || any(ynew <= 0) ) {
    dbt <- 0 #Outside simplex, always reject
  } else {
    dbt <- sum(dnorm(y, ynew, s, log = TRUE) - dnorm(ynew, y, s, log = TRUE))
  }
  #Return new logit-scaled point and corresponding dbt
  attr(ynew, 'dbt') <- dbt
  return(ynew)
}

## ----resPropSamp---------------------------------------------------------
RunMhRP <- function(center, B, concentration, s, type, dat = NULL){
  #Redefining parameters, setting initial time, and making empty vectors and matrices
  zz = proc.time();
  p <- length(center)
  S <- array(0, c(B, p)) 
  moveCount <- 0
  center <- center/sum(center)
  a <- concentration*center 
  ycurrent <- center
  #Run sampler
  for (i in 1:B) { 
    #Propose new y
    ycand <- PropStepRP(ycurrent, s)
    #Decide to accept or reject
    if(sum(ycand) != 1 || any(ycand <= 0)){
      S[i, ] <- ycurrent #reject if outside simplex
    } else {
      move <- (log(runif(1)) < attr(ycand, 'dbt') + TargetDir(ycand, ycurrent, a, dat))
      if (!is.na(move) & move){
        ycurrent <- ycand
        moveCount <- moveCount + 1
      } 
      #Store y for this iteration
      S[i, ] <- ycurrent
    }
  }
  #Timing
  runTime <- proc.time()-zz
  #Return results
  return(list(S = S, runTime = runTime, moveCount = moveCount, p = p, 
              center = center, B = B, concentration = concentration, s = s,
              type = type, dat = NULL, a = a))
}

## ----resPropTab----------------------------------------------------------
library(coda) #for effective sample size function
#Vector of possible s values
sVec <- c(0.0005, 0.001, 0.005, 0.01, 0.05, 0.1)
#Make empty table for results
tab2 <- matrix(nrow = 6, ncol = 3)
colnames(tab2) <- c("Step Size(s)", "Acceptance Rate", "Mean Effective Sample Size")
tab2[, 1] <- sVec
#Run sampler for each s value and record acceptance rate and mean effective sample size
for(i in 1:6){
  s <- rep(sVec[i], 20)
  highDim <- RunMhRP(center = rep(1/20, 20), B = 5e3, concentration = 20, s = s,
                    type = 'dirichlet', dat = NULL)
  tab2[i, 2] <- round(highDim$moveCount/highDim$B, 5)
  tab2[i, 3] <- round(mean(effectiveSize(highDim$S)), 5)
}
tab2

## ----thinning, fig.height = 4.5, fig.width = 4, fig.align = "center"-----
#Run chain
thinning <- RunMhRP(center = c(1/3, 1/3, 1/3), B = 1e4, concentration = 3,
                  s = rep(.3, 3), type = 'dirichlet', dat = NULL)
#Plot
TriPlot(thinning) #Use SALTSampler plotting function
mtext(sprintf("h = (2, 2, 2) Proposal, Acceptance Rate: %s", 
              round(thinning$moveCount/thinning$B, 3)[1], 
              round(thinning$moveCount/thinning$B, 3)[2],
              round(thinning$moveCount/thinning$B, 3)[3]), side = 1, line = 0, cex = 0.8)

## ----nothinnig, fig.height = 4.5, fig.width = 4, fig.align = "center"----
#Run chain
noThinning <- RunMh(center = c(1/3, 1/3, 1/3), B = 1e4, concentration = 3, h = c(2, 2, 2), 
                    type = 'dirichlet')
#Plotting
TriPlot(noThinning) 
mtext(sprintf("h = (2, 2, 2) Proposal, Acceptance Rates: (%s, %s, %s)",
              round(noThinning$moveCount/noThinning$B, 3)[1], 
              round(noThinning$moveCount/noThinning$B, 3)[2], 
              round(noThinning$moveCount/noThinning$B, 3)[3]), side = 1, line = 0, cex = 0.8)

## ----SALThighDim---------------------------------------------------------
#Run chain
highDim <- RunMh(center = rep(1/20, 20), B = 5e3, concentration = 20, h = rep(2.4, 20),
                 type = 'dirichlet')
#Mean acceptance rate
mean(highDim$moveCount/highDim$B)
#Mean effective sample size
mean(effectiveSize(highDim$Y)) 

## ----SALTordMag----------------------------------------------------------
#Run chain
ordMag <- RunMh(center = c(.0001, .01, .9899), B = 5e3, concentration = 1e6, 
                h = c(.2, .02, .02), type = 'dirichlet')
#Acceptance rates
ordMag$moveCount/ordMag$B
#Posterior means
apply(ordMag$S, 2, mean) 

## ----calibFn-------------------------------------------------------------
#Function to calibrate
CalibFn <- function(y, logit = FALSE) {
  if (logit == TRUE) {
    y <- exp(LogPq(y)$logp)
  }
  out <- 1e3*y[1]^3*y[2]^3/sqrt(20 + y[3])
  return(out)
}
#Generate samples
z <- rnorm(1000, CalibFn(c(1/3, 1/3, 1/3), 2))

## ----TargetFn------------------------------------------------------------
Target <- function(ycand, ycurrent, a, dat, pars = NULL) {
  out <- sum(dnorm(dat, CalibFn(ycand, logit = TRUE), 2, log = TRUE)) - 
         sum(dnorm(dat, CalibFn(ycurrent, logit = TRUE), 2, log = TRUE)) + 
         sum((a - 1)*(LogPq(ycand)$logp - LogPq(ycurrent)$logp))
  return(out)
} 

## ----inputDist-----------------------------------------------------------
inputDist <- RunMh(center = c(1/3, 1/3, 1/3), B = 3e4, concentration = 3, 
                   h = c(0.2, 0.2, 0.2), type = 'user', dat = z)

## ---- fig.height = 4.5, fig.width = 4, fig.align = "center"--------------
TriPlot(inputDist)
mtext("Projected Samples of Y", side = 3, cex = 1.25, font = 2)

## ---- fig.height = 4, fig.width = 7, fig.align = "center"----------------
par(mfrow = c(1,2))
Diagnostics(inputDist)

