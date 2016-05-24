pnew.pold <- function(vec, df, m, Mnew, Mold, cc, c.pnew, c.pold)
# This function is called by function lr() to compute R-N derivative
# (ratio of densities) for uncond. Dirichlet models with
# Dir. precision parameters Mnew (numerator) v. Mold (vector,
# denom.)  Normal or t models (numerator) v. normal model
# (denominators) with the same mu, tau parameters for numer. and
# denom. models (mu and tau are the centrality and dispersion
# parameters, resp.)
# This is for the left side of Eqn. (2.12), one MCMC simulation.
# The function lr() calls pnew.pold() repeatedly for each row of MCMC
# output, using the apply() function.
# In fact, we are using the improved, multi-chain algorithm of
# Doss (2009) Eqns. (2.5) and (2.6).
# vec looks like: rownumber, psi_1, ..., psi_m,
# mu, tau, int(XXX?)  df df parameter of t distr. center of Dir.
# model for numerator
  
  {
    k <- length(cc)
    psi <- vec[1:m]
    mu <- vec[m+1]
    tau <- vec[m+2]
    psi.uniq <- unique(psi)
    d <- length(psi.uniq)
    psi.zscore <- (psi.uniq - mu)/tau
    if (df == -99)
      { pnew.term1.vec <- dnorm(psi.zscore) }
    else
      { pnew.term1.vec <- dt(psi.zscore,df) }
    if (Mnew == -99)
      {
        if (d < m) pnew <- 0
        else  pnew <- prod(pnew.term1.vec)
      }
    else
      {
        pnew <- prod(pnew.term1.vec) * Mnew^d * c.pnew #a scalar
      }
    pold.term1.vec <- dnorm(psi.zscore)
    pold.term1 <- prod(pold.term1.vec)
    if (Mnew == -99)
      {
        if (d == m)  pold.vec <- pold.term1 * c.pold
        else pold.vec <- rep(1,k) #dummy for case Mnew == -99 & d < m
      }
    else
      {
        pold.vec <- pold.term1 * Mold^d * c.pold
      }
    pold <- (t(cc) %*% pold.vec)/k
    lr <- pnew/pold
    lr
  }

lr <- function(df=-99, Mnew, Mold, cc, mat.list)
# computes the Bayes factors for uncond. Dir. models with the same psi, tau
# one Mnew (num.) and a vector of Mold (den.)
# Input:
#     mat.list a list containing k matrices of MCMC output
#              from uncond. Dir. model
# 
  {
    k <- length(cc) #no. of values of Mo
    nc <- ncol(mat.list[[1]]) #The matrices have same no. cols
    m <- nc-3
    if (Mnew == -99)
      {
        c.pnew <- 1 #dummy value
        log.c.pold <- lgamma(Mold + 1) - lgamma(Mold + m)
        c.pold <- Mold ^ (m-1) * exp(log.c.pold)
      }
    else
      {
        log.c.pnew <- lgamma(Mnew) - lgamma(Mnew+m)
        c.pnew <- exp(log.c.pnew)
        log.c.pold <- lgamma(Mold) - lgamma(Mold + m)
        c.pold <- exp(log.c.pold)
      }
    answer.list <- list()
    for (i in 1:k)
      {
       answer.list[[i]] <- apply(mat.list[[i]],1,pnew.pold,
                                 df=df,m=m,Mnew=Mnew,Mold=Mold,
                                 cc=cc,c.pnew=c.pnew,c.pold=c.pold)
     }
    unlist(answer.list)
  }


pnew.pold.c.o <- function(vec, df, m, Mnew, Mold, cc,log.c.pnew, c.pold)
  {
    k <- length(cc)
    psi <- vec[1:m]
    mu <- vec[m+1]
    tau <- vec[m+2]
    psi.uniq <- unique(psi)
    d <- length(psi.uniq)
    mminus <- sum(psi <= mu)
    mplus <- m - mminus
    psi.zscore <- (psi.uniq - mu)/tau
    if (df == -99)
      { pnew.term1.vec <- dnorm(psi.zscore) }
    else
      { pnew.term1.vec <- dt(psi.zscore,df) }
    if (Mnew == -99)
      {
        if (d < m) pnew <- 0
        else  pnew <- prod(pnew.term1.vec)
      }
    else
      {
        Mn2 <- Mnew/2
        log.c2 <- log.c.pnew - lgamma(Mn2 + mminus) - lgamma(Mn2 + mplus) 
        pnew <- prod(pnew.term1.vec) * Mnew^d * exp(log.c2)  #a scalar
      }
    pold.term1.vec <- dnorm(psi.zscore)
    pold.term1 <- prod(pold.term1.vec)
    if (Mnew == -99)
      {
        if (d == m)  pold.vec <- pold.term1 * c.pold
        else pold.vec <- rep(1,k) #dummy for case Mnew == -99 & d < m
      }
    else
      {
        pold.vec <- pold.term1 * Mold^d * c.pold
      }  
    pold <- (t(cc) %*% pold.vec)/k
    lr <- pnew/pold
    lr
  }

lr.c.o <- function(df=-99, Mnew, Mold, cc, mat.list)
  # mat.list  list of matrices of MCMC output from uncond. Dir.
  #           See bottom of p. 13
  {
    k <- length(cc)
    nc <- ncol(mat.list[[1]]) #The matrices have same no. cols
    m <- nc-3
    if (Mnew == -99)
      {
        c.pnew <- 1 #dummy value
        log.c.pold <- lgamma(Mold + 1) - lgamma(Mold + m)
        c.pold <- Mold ^ (m-1) * exp(log.c.pold)
      }  
    else
      {  
        log.c.pnew <- -m *log(2) + 2*lgamma(Mnew/2)
#        c.pnew <- 1/2^m * exp(log.c.pnew)
        log.c.pold <- lgamma(Mold) - lgamma(Mold + m)
        c.pold <- exp(log.c.pold)
      }
    answer.list <- list()
    for (i in 1:k)
      {
       answer.list[[i]] <- apply(mat.list[[i]],1,pnew.pold.c.o,
                                 df=df,m=m,Mnew=Mnew,Mold=Mold,
                                 cc=cc,log.c.pnew=log.c.pnew,c.pold=c.pold)
     }
    unlist(answer.list)
  }    


# m1.c/m div. by m1.o/m  =  m1.c/m1.o, the Bayes factor for cond. v.
# uncond. Dir. models with Dir. precision parameter Mnew
bf.c.o <- function(df=-99, from=.4, incr=.1, to, cc, mat.list)
# function to create object to plot; this is slow
  {
   if (missing(to)) stop("argument 'to' is missing")
   if (to < .4 || to > 100) stop("'to' must be > .4 and < 100")
   if (!is.numeric(cc) || length(cc) != 9){
      stop("'cc' must be numeric vector length 9")
    }
    nM <- length(mat.list)
    if (!is.list(mat.list) && nM != 9) {stop("'mat.list' must be a list
      of 9 matrices of MCMC output")}
    dOut <- dim(mat.list[[1]])
    ncycles <- dOut[1] - 1
    nparam <- dOut[2]
    for (ii in 1:nM) 
      {
        dOutii <- dim(mat.list[[ii]])
        if (dOutii[1] -1 != ncycles || dOutii[2] != nparam)
          {stop("all matrices in mat.list must have same size\n")}
      }    
    Mold <- c(.25,.5,1,2,4,8,16,32,64)
    Mnew <- c(seq(from=from, to=to/5, by=incr),
              seq(from=to/5, to=to, by=5*incr)[-1])
    y <- vector(length=length(Mnew))
    for (i in 1:length(Mnew))
    {
      y[i] <-
        mean(lr.c.o(df=df, cc=cc, Mnew=Mnew[i], Mold=Mold, mat.list=mat.list)) /
        mean(lr(df=df, cc=cc, Mnew=Mnew[i], Mold=Mold, mat.list=mat.list))
    }
    list(Mnew=Mnew,y=y,yinfinity=1)
  }

draw.bf <- function(obj,line.lwd=2, ...)
{
  if (missing(obj)) stop("object to plot is missing")
  ncomp <- length(obj)
  if (!is.list(obj) || length(obj) != 3) stop("obj needs 3 components")
  Mnew <- obj[[1]]
  y <- obj[[2]]
  yi <- obj[[3]]
  if (length(Mnew) != length(y) || length(yi) > 1) stop("lengths
 of the two vectors must match, and yi must be a single value")
  oldpar <- par(no.readonly=TRUE)  
  ylim.u <- max(y) + .5
  to <- max(Mnew)/.88
  oldpar.maretc <- par(mar=c(4,4,0.5,.5),...)#c(bottom, left, top, right)
# default is 'c(5, 4, 4, 2) + 0.1'.
#  par(cex=.75)
  plot(0.1,0.1, type="n", xlab="", ylab="",
       xlim=c(0,to),  ylim=c(0,ylim.u), xaxt="n", bty="n",...)
  axis(1,at=seq(0, to*.88, by=1))
# par(cex=1.4); axis(1, at=to, expression(infinity), tck=0); par(cex=.75)
  oldpar.cex <- par(cex=1.4)
  axis(1, at=to, expression(infinity), tck=0)
  par(cex=oldpar.cex)
#  par(cex=1); title(ylab="Bayes Factor")
  title(ylab="Bayes Factor")
#  par(cex=1); title(xlab="M")
  title(xlab="M")
  from <- .4
#  par(lwd=4)
#  oldpar.lwd <- par("lwd")
  lines(Mnew, y, type="l", lwd=line.lwd)
#  lines(Mnew, y, type="l")
#  par(lwd=4-1)
  points(to, yi, lwd=max(1,line.lwd-1))
  par(oldpar)
}


bf.o <- function(df=-99, from=.4, incr=.1, to, cc, mat.list)
# function to create object to plot; this is slow
  {
   if (missing(to)) stop("argument 'to' is missing")
   if (to < .4 || to > 100) stop("'to' must be > .4 and < 100")
   if (!is.numeric(cc) || length(cc) != 9){
      stop("'cc' must be numeric vector length 9")
    }
    nM <- length(mat.list)
    if (!is.list(mat.list) || nM != 9) {stop("'mat.list' must be a list
      of 9 matrices of MCMC output")}
    dOut <- dim(mat.list[[1]])
    ncycles <- dOut[1] - 1
    nparam <- dOut[2]
    for (ii in 1:nM) 
      {
        dOutii <- dim(mat.list[[ii]])
        if (dOutii[1] -1 != ncycles || dOutii[2] != nparam)
          {stop("all matrices in mat.list must have same size\n")}
      }    
    Mold <- c(.25,.5,1,2,4,8,16,32,64)
    Mnew <- c(seq(from=from, to=to/5, by=incr), seq(from=to/5, to=to, by=5*incr))
    lM <- length(Mnew)
    y <- vector(length=lM)
    for (i in 1:lM)
    {
      y[i] <- mean(lr(df=df, Mnew=Mnew[i],cc=cc,
                      Mold=Mold, mat.list=mat.list))
    }
  yinfinity <- mean(lr(Mold=Mold, df=-99,Mnew=-99,cc=cc,mat.list=mat.list))  
  list(Mnew=Mnew,y=y,yinfinity=yinfinity)
  }

bf.c <- function(df=-99, from=.4, incr=.1, to, cc, mat.list)
# function to create object to plot; this is slow
  {
   if (missing(to)) stop("argument 'to' is missing")
   if (to < .4 || to > 100) stop("'to' must be > .4 and < 100")
   if (!is.numeric(cc) || length(cc) != 9){
      stop("'cc' must be numeric vector length 9")
    }
    nM <- length(mat.list)
    if (!is.list(mat.list) && nM != 9) {stop("'mat.list' must be a list
      of 9 matrices of MCMC output")}
    dOut <- dim(mat.list[[1]])
    ncycles <- dOut[1] - 1
    nparam <- dOut[2]
    for (ii in 1:nM) 
      {
        dOutii <- dim(mat.list[[ii]])
        if (dOutii[1] -1 != ncycles || dOutii[2] != nparam)
          {stop("all matrices in mat.list must have same size\n")}
      }    
    Mold <- c(.25,.5,1,2,4,8,16,32,64)
    Mnew <- c(seq(from=from, to=to/5, by=incr), seq(from=to/5, to=to, by=5*incr))
    lM <- length(Mnew)
    y <- vector(length=lM)
    for (i in 1:lM)
    {
      y[i] <- mean(lr.c.o(df=df, Mnew=Mnew[i],cc=cc, Mold=Mold, mat.list=mat.list))
    }
  yinfinity <- mean(lr.c.o(Mold=Mold, df=-99,Mnew=-99,cc=cc,mat.list=mat.list))  
  list(Mnew=Mnew,y=y,yinfinity=yinfinity)
  }
