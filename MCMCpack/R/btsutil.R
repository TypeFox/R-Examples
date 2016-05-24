##########################################################################
## Utility Functions for Bayesian Times Series Models
##
## written and maintained by:
##    Jong Hee Park
##    Department of Political Science
##    University of Chicago
##    jhp@uchicago.edu
##
## Revised on 09/12/2007 JHP	  
##
## NOTE: only the plot functions are documented and exported in the
## NAMESPACE.
##
## This software is distributed under the terms of the GNU GENERAL
## PUBLIC LICENSE Version 2, June 1991.  See the package LICENSE
## file for more information.
##
## Copyright (C) 2003-2007 Andrew D. Martin and Kevin M. Quinn
## Copyright (C) 2007-present Andrew D. Martin, Kevin M. Quinn,
##    and Jong Hee Park
##########################################################################



##############################################################
## Helper functions for MCMCpoissonChange and MCMCbinaryChange()
##############################################################
    
## switch a state vector into a matrix containing the number of states
"switchg" <-  function(s1){    
  s <- max(s1)
  out <- matrix(0,s,s)
  
  ## all P(i,i+1) are 1
  for (i in 1:(s-1)){
    out[i,i+1] <- 1}
  
  ## diagonal elements is (the number of occurrence - 1)
  diag(out) <- table(s1)-1
  return(out)
}

## "trans.mat.prior" makes a transition matrix
"trans.mat.prior" <- function(m, n, a=NULL, b=NULL){
  if (!is.null(a)|!is.null(b)){
    a <- a
    b <- b
  }
  else {
    expected.duration <- round(n/(m+1))
    b <- 0.1
    a <- b*expected.duration
  }
  trans <- diag(1, m+1)
  ## put a as diagonal elements except the last row
  diag(trans)[1:m]<-rep(a, m)
  ## put b in trans[i, i+1]
  for (i in 1:m){trans[i, i+1]<-b}
  return(trans)
} 

## "plotState" draws a plot of posterior distribution of states 
"plotState" <-
  function (mcmcout, main="Posterior Regime Probability", ylab=expression(paste("Pr(", S[t], "= k |", Y[t], ")")),
            legend.control = NULL, cex = 0.8, lwd = 1.2, start=1)
  {
    out <- attr(mcmcout, "prob.state")
    y <- attr(mcmcout, "y")
    m <- attr(mcmcout, "m")
    
    if (!is.ts(y))
      y <- ts(y, start)
    time.frame <- as.vector(time(y))
    
    plot(start, 0, xlim = range(time.frame), ylim = c(0, 1), type = "n",
         main = main, xlab = "Time", cex = cex, lwd = lwd,
         ylab = ylab, axes=F)
    axis(1, tick = FALSE, col="darkgrey")
    axis(2, tick = FALSE, col="darkgrey")
    for (i in 1:length(axTicks(1))) lines(c(axTicks(1)[i], axTicks(1)[i]), c(0,1), col="darkgrey")
    for (i in 1:length(axTicks(2))) lines(c(axTicks(2)[i], max(axTicks(1))), c(axTicks(2)[i], axTicks(2)[i]), col="darkgrey")
     
    for (i in 1:(m + 1)) points(time.frame, out[, i], type = "o",
                                lty = i, lwd = lwd, col = i, cex = cex)
    if (!is.null(legend.control)) {
      if (length(legend.control) != 2)
        stop("You should specify x and y coordinate for a legend.")
      else legend(legend.control[1], legend.control[2],
                  legend = paste("State",1:(m + 1), sep = ""),
                  col = 1:(m + 1), lty = 1:(m + 1),
                  lwd = rep(lwd, m + 1), pch = rep(1, m + 1), bty = "n")
    }
  }


## "plotChangepoint" draws a plot of posterior changepoint probability
## Thanks to Patrick Brandt for providing the idea of overlaying.
"plotChangepoint" <-
  function (mcmcout, main="Posterior Density of Regime Change Probabilities", xlab = "Time", ylab = "",
            verbose = FALSE, start=1, overlay=FALSE)
  {
    out <- attr(mcmcout, "prob.state")
    y <- attr(mcmcout, "y")
    m <- attr(mcmcout, "m")
    if(overlay==FALSE){
      par(mfrow = c(m, 1), mar = c(2, 4, 1, 1))
    }
    if (!is.ts(y))
      y <- ts(y, start)
    time.frame <- as.vector(time(y))
    
    if (m == 1) {
      pr.st <- c(0, diff(out[, (m + 1)]))
      pr.st[pr.st<0] <- 0           
      plot(time.frame, pr.st, type = "h", lwd=2, main = main, xlab = xlab, ylab = ylab, axes=F)
      axis(1, tick = FALSE, col="darkgrey")
      axis(2, tick = FALSE, col="darkgrey")
      for (i in 1:length(axTicks(1))) lines(c(axTicks(1)[i], axTicks(1)[i]), c(0, max(axTicks(2))), col="darkgrey")
      for (i in 1:length(axTicks(2))) lines(c(axTicks(2)[i], max(axTicks(1))), c(axTicks(2)[i], axTicks(2)[i]), col="darkgrey")
      cp <- which(cumsum(pr.st) > 0.5)[1] - 1
      lines(c(cp + time.frame[1], cp + time.frame[1]), c(0, max(axTicks(2))), lty = 3, col = "red")
    }
    
    else {
      cp <- rep(NA, m)
      for (i in 2:m) {
        pr.st <- c(0, diff(out[, i]))
        pr.st <- ifelse(pr.st < 0, 0, pr.st)
        plot(time.frame, pr.st, type = "h", lwd=2, main = "", xlab = xlab, ylab = ylab, col="black", axes=FALSE)
        axis(1, tick = FALSE, col="darkgrey")
        axis(2, tick = FALSE, col="darkgrey")
        for (k in 1:length(axTicks(1))) {lines(c(axTicks(1)[k], axTicks(1)[k]), c(0, max(axTicks(2))), col="darkgrey")}
        for (k in 1:length(axTicks(2))) {lines(c(axTicks(2)[k], max(axTicks(1))), c(axTicks(2)[k], axTicks(2)[k]),
                                              col="darkgrey")}
        cp[i - 1] <- which(cumsum(pr.st) > 0.5)[1] - 1
        lines(c(cp[i - 1] + time.frame[1], cp[i - 1] + time.frame[1]), c(0, max(axTicks(2))), lty = 3, col = "red")
      }
      pr.st <- c(0, diff(out[, (m + 1)]))
      pr.st[pr.st<0] <- 0           
      plot(time.frame, pr.st, type = "h", lwd=2, main = main, xlab = xlab, ylab = ylab, col="black", axes=FALSE)
      axis(1, tick = FALSE, col="darkgrey")
      axis(2, tick = FALSE, col="darkgrey")
      for (k in 1:length(axTicks(1))) {lines(c(axTicks(1)[k], axTicks(1)[k]), c(0, max(axTicks(2))), col="darkgrey")}
      for (k in 1:length(axTicks(2))) {lines(c(axTicks(2)[k], max(axTicks(1))), c(axTicks(2)[k], axTicks(2)[k]),
                                              col="darkgrey")}
      cp[m] <- which(cumsum(pr.st) > 0.5)[1] - 1
      lines(c(cp[m] + time.frame[1], cp[m] + time.frame[1]), c(0, max(axTicks(2))), lty = 3, col = "red")
    }

    cp.means <- rep(NA, m + 1)
    cp.start <- c(1, cp + 1)
    cp.end <- c(cp, length(y))

    if (verbose == TRUE){
      cat("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n")
      cat("Expected changepoint(s) ", cp + time.frame[1], "\n")
      for (i in 1:(m + 1)) cp.means[i] <- mean(y[cp.start[i]:cp.end[i]])
      cat("Local means for each regime are ", cp.means, "\n")
      cat("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n")
    }
  }

## prior check for transition matrix
"check.P" <- function(P.start = NA, m=m, n=n, a=a, b=b){
  if (is.na(P.start)[1]){
    P <- trans.mat.prior(m=m, n=n, a=0.9, b=0.1)}
  else if ((dim(P.start)[1]==m+1)&(dim(P.start)[2]==m+1)){
    if ((max(P.start)>1)||(min(P.start)<0)){
      stop("Error: P starting values are out of unit range.\n")
    }
    else
      P <- P.start
  }
  else {
    stop("Error: P starting values are not conformable.\n")
  }
  return(P)
}

## priro check for mean
"check.theta" <- function(theta.start = NA, ns = ns, y = y, min, max){
  if (is.na(theta.start)[1]){ # draw from uniform with range(y)
    theta <- runif(ns, min=min, max=max)}
  else if (length(theta.start)==ns)
    theta <- theta.start
  else if (length(theta.start)!=ns) {
    stop("Error: theta starting values are not conformable.\n")
  }
  return(theta)
}


## initial values of tau in MCMCpoissonChange
"tau.initial" <- function(y, tot.comp){
  tau             <-  rep(NA, tot.comp)
  lambda.t        <-  0.1
  count           <-  0
  for (t in 1:length(y)){
    nt      <-  y[t]
    if (nt==0) {
      taut    <-  1 + rexp(1, lambda.t)
        count   <-  count + nt + 1
    }
    else{
      ut      <-  runif(nt)
      uorder  <-  c(0, sort(ut))
      tau.tj  <-  diff(uorder)
      sum.tau.tj <- sum(tau.tj)
      tau.last<-  1 - sum.tau.tj + rexp(1, y[t])
      count   <-  count + nt + 1
      taut    <-  c(tau.tj, tau.last)
    }
    tau[(count-nt):count] <-   taut
  }
  return(tau)
}


## beta starting values in MCMCpoissonChange()
"beta.change.start"<- function (beta.start, ns, k, formula, family, data){
  ## if a user does not specify beta.start, use a coefficient vector from mle
  if (is.na(beta.start[1])) {
    b0 <- coef(glm(formula, family = family,  data = data))
    beta.start  <-  matrix(rep(b0, ns), ns, k, byrow=TRUE)
  }
  ## if beta.start is scalar or k by 1 vector, repeat this
  else if (is.null(dim(beta.start))&&length(beta.start)<=k) {
    beta.start <- beta.start * matrix(1, ns, k)
    ## this alternates beta.start if beta.start is not a scalar
  }
  ## if the length of beta.start is same to ns*k, make this as a matrix
  else if (is.null(dim(beta.start))&&length(beta.start)==ns*k) {
    beta.start <- matrix(beta.start, ns, k)
  }
  else if (is.null(dim(beta.start))&&length(beta.start)>=k) {
    cat("Error: Starting value for beta not conformable.\n")
    stop("Please respecify and call ", calling.function(),
         " again.\n", call. = FALSE)
  }
  ## else, report an error message and stop
  else if (!all(dim(beta.start) == c(ns, k))) {
    cat("Error: Starting value for beta not conformable.\n")
    stop("Please respecify and call ", calling.function(),
         " again.\n", call. = FALSE)
  }
  return(beta.start)
}


## draw predicted outcomes of intervention analysis
plotIntervention <- function(mcmcout, forward = TRUE, start = 1,
                             alpha = 0.05, col = "red", main="", ylab="", xlab=""){
  ## pull inputs
  y <- ts(attr(mcmcout, "y"), start=start)
  inter <- attr(mcmcout, "intervention") 
  N <- length(y)
  
  ## draw a plot
  plot(y, main=main, ylab=ylab, xlab=xlab)
  abline(v = time(y)[inter], lty=2)
  if (forward == TRUE){
    yfore <- attr(mcmcout, "yforepred")
    yfore.mu <- ts(apply(yfore, 2, mean), start=start); yfore.mu[1:(inter-1)] <- NA
    yfore.upper <- ts(apply(yfore, 2, quantile, probs=(1-alpha/2)), start=start); yfore.upper[1:(inter-1)] <- NA
    yfore.lower <- ts(apply(yfore, 2, quantile, probs=(alpha/2)), start=start); yfore.lower[1:(inter-1)] <- NA
    lines(yfore.mu, col=col, lwd=2)
    lines(yfore.upper, col=col, lty=3)
    lines(yfore.lower, col=col, lty=3)
  }
  else {
    yback <- attr(mcmcout, "ybackpred")
    yback.mu <- ts(apply(yback, 2, mean), start=start);  yback.mu[(inter+1):N] <- NA
    yback.upper <- ts(apply(yback, 2, quantile, probs=(1-alpha/2)), start=start); yback.upper[(inter+1):N] <- NA
    yback.lower <- ts(apply(yback, 2, quantile, probs=(alpha/2)), start=start); yback.lower[(inter+1):N] <- NA
    lines(yback.mu, col=col, lwd=2)
    lines(yback.upper, col=col, lty=3)
    lines(yback.lower, col=col, lty=3)
  }
}
## Example
## pdf(file="Nile_MCMCinter.pdf", width=12, height=4)
## par(mfrow=c(1,3))
## plotState(ar1, start=1871, main="Hidden Regime Change")
## plotIntervention(ar1, start=1871, main="Forward Analysis", alpha= 0.5, ylab="Nile River flow", xlab="Year")
## plotIntervention(ar1, forward=FALSE, start=1871, main="Backward Analysis", alpha= 0.5, ylab="Nile River flow", xlab="Year")
## dev.off()

## when we compare models in a list
BayesFactorList <- function (model.list){
  oldM <- length(model.list)
  zero.marg <- rep(NA, oldM)
  for (j in 1:oldM) {
    zero.marg[j] <- ifelse(attr(model.list[[j]], "logmarglike") == 0, 1, 0)
  }
  new.model.list <- model.list[c(which(zero.marg == 0))]
  M <- length(new.model.list)
  out <- matrix(NA, M, 2)
  BF <- rep(NA, M)
  for (j in 1:M) {
    BF[j] <- attr(new.model.list[[j]], "logmarglike")
    out[j, 1] <- BF[j]
  }
  if (sum(exp(BF) == 0)){
    ## if log like is too small, add some constants
    BF <- BF + abs(BF[1])
  }
  prob <- exp(BF)/sum(exp(BF))
  out[, 2] <- prob
  marker <- which(zero.marg == 0)
  rownames(out) <- names(model.list)[marker]
  return(out)
}

## consecutive geweke diag test
geweke.test <- function(output.list, z.value=1.96){
  n <- length(output.list)
  result <- rep(NA, n)
  cat("\n --------------------------------- ")
  for (i in 1:n){
    if(sum(abs(geweke.diag(output.list[[i]])$z) > z.value)>0){
      cat("\n Non-convergence for model ", i)
      result[i] <- "Fail"
    }
    else {
      result[i] <- "Pass"
    }
  }
  cat("\n --------------------------------- \n")
  return(result)
}
## outputlist <- list(ar0, ar1a, ar2a, ar1f, ar2f, ar1r, ar2r, tr0, tr1a, tr2a, tr1f, tr2f, tr1r, tr2r)
## conv <- geweke.test(outputlist)

## consecutive heidel Heidelberger and Welch's convergence diagnostic test
heidel.test <- function(output.list, p.value=0.05){
  n <- length(output.list)
  result <- rep(NA, n)
  cat("\n --------------------------------- ")
  for (i in 1:n){
    print(i)
    plist1 <- heidel.diag(output.list[[i]], pvalue=p.value)[,1]
    plist2 <- heidel.diag(output.list[[i]], pvalue=p.value)[,4]
    if(sum(c(plist1, plist2) == 0)>0){
      cat("\n Non-convergence for model ", i)
      result[i] <- "Fail"
    }
    else {
      result[i] <- "Pass"
    }
  }
  cat("\n --------------------------------- \n")
  return(result)
}
## outputlist <- list(ar0, ar1a, ar2a, ar1f, ar2f, ar1r, ar2r, tr0, tr1a, tr2a, tr1f, tr2f, tr1r, tr2r)
## conv <- geweke.test(outputlist)
