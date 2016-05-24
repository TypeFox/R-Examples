# $Id: dpqperm.R,v 1.24.4.1 2004/09/27 06:07:11 hothorn Exp $

findfact <- function(scores, m, tol) {
    sc <- rev(sort(scores))
    foo <- function(fact) {
      abs(sum(abs(fact*sc[1:m] - round(fact*sc[1:m],0)))/fact - tol)
    }
    grid <- seq(from = 1, to = 10000/(sum(sc) -min(sc)), by = 1)
    fact <- which.min(sapply(grid, foo))
    ## fact <- optimize(foo, c(1, 100000/(sum(sc) -min(sc))))$minimum
    thistol <- sum(abs(fact*sc[1:m] - round(fact*sc[1:m],0)))/fact
    if (thistol > tol)
      warning(paste("cannot hold tol, tolerance:", round(thistol, 6)))
    fact
}


dperm <- function(x, scores, m, paired = NULL, tol = 0.01, fact = NULL,
                  density=FALSE, simulate=FALSE, B=10000)
{
    if (is.null(x)) stop("Non-numeric argument to mathematical function")
    # use raw scores for Monte-Carle simulations
    eq <- equiscores(scores, m, tol, fact, simulate)
    if (simulate) 
      cp <- cperm(eq, m, paired, B=B)
    else 
      cp <- cperm(eq, m, paired)
    # return everything up to x
    if (density & length(x) == 1) {
      indx <- which(cp$T <= x)
      RVAL <- data.frame(cp$T[indx], cp$Prob[indx])
    } else {
      RVAL <- rep(0, length(x))
      RVAL[x %in% cp$T ] <- cp$Prob[cp$T %in% x]
    }
    return(RVAL) 
}


pperm <- function(q, scores, m, paired = NULL, tol = 0.01, fact = NULL,
                     alternative=c("less", "greater", "two.sided"),
                     pprob=FALSE, simulate=FALSE, B=10000)
{
    if(is.null(q)) stop("Non-numeric argument to mathematical function")
    alternative <- match.arg(alternative)
    if (is.null(paired))  paired <- (length(scores) == m)
    eq <- equiscores(scores, m, tol, fact, simulate)
    if (simulate) 
      cp <- cperm(eq, m, paired, B=B)
    else 
      cp <- cperm(eq, m, paired)
    PVALUE <- c()
    PPROB <- c()
    for (i in q) {
        if(alternative=="less")
            prob <- cp$Prob[cp$T <= i]
        if(alternative=="greater")
            prob <- cp$Prob[cp$T >= i]
        if(alternative=="two.sided" & paired) {
            if(2*sum(cp$Prob[cp$T <= i]) <= 1)
                prob <- 2*cp$Prob[cp$T <= i]
            else
                prob <- 2*cp$Prob[cp$T >= i]
        }
        if(alternative=="two.sided" & !paired) {
            expect <- m/length(scores)*sum(scores)
            TH <- cp$T - expect
            ih <- i - expect
            if (ih != 0) 
              prob <- c(cp$Prob[TH <= ifelse(ih > 0, -ih, ih)],
                        cp$Prob[TH >= ifelse(ih >= 0, ih, -ih)])
            else
              prob <- 1
        }
    PVALUE <- c(PVALUE, min(c(1, sum(prob)))) ### STATISTIC == Expectation
    PPROB <- c(PPROB, ifelse(!is.null(cp$Prob[cp$T == i]),
                                      cp$Prob[cp$T ==i], 0)) 
    }
    if(pprob)
        return(list(PVALUE = PVALUE, PPROB = PPROB))
    else
        return(PVALUE)
} 


qperm <- function(p, scores, m, paired = NULL, tol = 0.01, fact = NULL,
                  simulate = FALSE, B=10000)
{
    if (is.null(p)) stop("Non-numeric argument to mathematical function")
    if (any(p < 0) || any(p > 1)) {
        warning("p is not a probability")
        return(NaN)
    }
    eq <- equiscores(scores, m, tol, fact, simulate)
    if (simulate) 
      cp <- cperm(eq, m, paired, B=B)
    else 
      cp <- cperm(eq, m, paired)
    cs <- cumsum(cp$Prob)
    RVAL <- c()
    for (i in p) {
        quant <- which(cs < i)
        if (length(quant) == 0) quant <- 0
        RVAL <- c(RVAL, cp$T[max(quant) + 1])
    }
    return(RVAL)
}

rperm <- function(n, scores, m)
    sapply(1:n, dummy <- function(x) sum(sample(scores,m)))

equiscores <- function(scores, m=NULL, tol = 0.01, fact=NULL, simulate=FALSE)
{
  if (any(is.null(scores))) 
      stop("Non-numeric argument to mathematical function")

  if (!is.null(m)) { 
    if (m < 1) 
        stop("m less than 1")
    if (m > length(scores)) 
        stop("m greater length(scores)")
    paired <- (length(scores) == m)
  } else {
    paired <- FALSE
  }

  # use raw scores for simulation based computation of the permutation 
  # distribution
  if (simulate) {
    if(paired & any(scores < 0)) stop("scores must be positive")
    RVAL <- list(scores = scores, fact = 1, add = 0)
    class(RVAL) <- "equis"
    return(RVAL)
  }  

  if (is.null(fact)) {
    # first, handle integer and midranked scores, 
    # m is not needed here
    fact <- 1
    fscore <- scores - floor(scores)
    if (any(fscore != 0)) {
      if (all(fscore[fscore != 0] == 0.5)) # midranks
        fact <- 2
    }
    scores <- scores*fact

    fscore <- scores - floor(scores)

    # ok, lets face the problems    
    if (!all(fscore == 0)) { # rational or real scores
      fact <- ceiling(findfact(scores - min(scores) + 1, m, tol))
      scores <- round(scores * fact)
    }
  } else {
    # the user knows what to do
    scores <- round(scores * fact)
  } 
   
  if(!paired) 
    add <- min(scores) - 1
  else 
    add <- 0
  scores <- scores - add

  if(any(scores < 0)) stop("scores must be positive")

  RVAL <- list(scores = scores, fact = fact, add = add)
  class(RVAL) <- "equis"
  return(RVAL)
}


cperm <- function(escores, m, paired = NULL, B=NULL)
{
    if (!(class(escores) == "equis"))
        stop("scores are not of class equis") 

    N <- length(escores$scores)

    if (is.null(paired))
        paired <- (N == m)
    else 
        paired <- (N == m) && paired  

    if (!is.null(B)) {
        if (B < 1) stop("B must be positive integer")
        RVAL <- .Call("sim2is", scores=as.double(escores$scores),
                                mfirst=as.integer(m), Nsim=as.integer(B), 
                                PACKAGE="exactRankTests")
        names(RVAL) <- c("T", "Prob")
        if (escores$add != 0 & paired) 
          warning("escores$add not zero for paired samples!")
        RVAL$T <- (RVAL$T + escores$add*m)/escores$fact
        return(RVAL)
    }

    if (paired) {
        # paired two sample situation
        prob <- .Call("cpermdist1", scores = as.integer(escores$scores),
                   PACKAGE="exactRankTests")
        t <- which(prob != 0)
        prob <- prob[t]
        # 0 is possible
	t <- t - 1
        if (escores$add != 0) 
          warning("escores$add not zero for paired samples!")
    } else {
        # independent samples
        col <- sum(sort(escores$scores)[(N + 1 - m):N])
        scores <- rep(1, N)
        prob <- .Call("cpermdist2", ma = as.integer(m),
                mb = as.integer(col), scorea = as.integer(scores), 
                scoreb = as.integer(escores$scores),
                retProb = as.logical(TRUE), PACKAGE="exactRankTests")
        t <- which(prob != 0)
        prob <- prob[t]
    }
    t <- (t + escores$add*m)/escores$fact
    RVAL <- list(T = t, Prob = prob)
    class(RVAL) <- "cperm"
    return(RVAL)
}
    