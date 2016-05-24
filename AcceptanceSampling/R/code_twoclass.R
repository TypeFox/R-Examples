## code_twoclass.R --- 
##
## Author: Andreas Kiermeier
##
## Created: 08 Mar 2007 
##
## Purpose: A package to provide functionality for creating and
##          evaluating acceptance sampling plans.
##          
## Changes:
## 16Aug07: * Added check in OC2c validation code to ensure that sample sizes
##            are greater than zero.
##          * Added virtual class OCvar for variables sampling plans (single
##            only)
##          * Added actual class for variables sampling plans - Normal
## 20Aug07: * Added function {find.k} to find constant k for given sample size in
##            normal variables sampling plans
##          * Added function {find.plan} to find smallest sampling plan
##            for given Producer and Consumer Risk Points
## 27Feb08: * Changed the validation for r & c (which are cumulative) to be
##            compared against cumsum(n) as n is not cumulative.
## 05Mar08: * Fixed problem with multiple sampling plans - previous calculations
##            were completely wrong.
##            Code now enumerates over all stages and all possible outcomes which
##            leadto additional sampling.
##            
## Notes:
## For implemented package use
## 
## getFromNamespace(paste("calc.",OCtype,sep=""), ns="AcceptanceSampling")
##
## while for testing directly use
##
## get(paste("calc.",OCtype,sep=""))
##
## There are THREE (3) of these instances
## ----------------------------------------------------------------------

## ----------------------------------------------------------------------
## Class definitions
## ----------------------------------------------------------------------

setClass("OC2c", representation(n="numeric", ## A vector of sample sizes at each
                                ## stage of sampling
                                ## NOT CUMULATIVE
                                c="numeric", ## vector of acceptance numbers for
                                ## each stage of sampling. Accept if actual number
                                ## of defectives/defects is <= c
                                ## CUMULATIVE
                                r="numeric", ## vector of rejection numbers for
                                ## each stage of sampling. Reject if actual number
                                ## of defectives/defects is >= r
                                ## CUMULATIVE
                              type="character",
                              paccept="numeric",
                              "VIRTUAL"),
         validity=function(object){
           if(any(is.na(object@n)) | any(is.na(object@c)) |
              any(is.na(object@r)))
             return("Missing values in 'n', 'c', or 'r' not allowed")
           ## Check that n, c and r are of the same length
           l <- length(object@n)
           if (l != length(object@c) | l != length(object@r))
             return("n, c and r must be of same length.")
           ## Check that the sample sizes make sense
           if (any(object@n <= 0))
             return("Sample size(s) 'n' must be greater than 0.")
           ## Check that the acceptance numbers make sense
           if (any(object@c < 0) | any(object@c > cumsum(object@n)))
             return("Acceptance number(s) 'c' must be in the range [0,n], here n is the cumulative sample size.")
           ## Check that the rejection numbers make sense
           if (any(object@r < 0) | any(object@r > cumsum(object@n)))
             return("Rejection number(s) 'r' must be in the range [0,n], here n is the cumulative sample size.")
           if (any(object@r <= object@c))
             return("Rejection number(s) 'r' must be greater than acceptance number(s) 'c'.")
           ## For double sampling (or more) make sure that acceptance and
           ## rejection number are non-decreasing and non-increasing, respectively.
           if (l > 1) {
             if (any(diff(object@c)<0) )
               return("'c' must be non-decreasing")
             if (any(diff(object@r)<0) )
               return("'r' must be non-decreasing")
           }
           ## Check that a decision is made on the last sample
           if (object@r[l] != object@c[l] + 1)
             return("Decision from last sample cannot be made: r != c+1")
           ## Otherwise things seem fine.
           return(TRUE)
         })

setClass("OCbinomial",
         representation("OC2c",
                        pd="numeric"),
         contains="OC2c",
         prototype=list("OC2c", type="binomial", pd=seq(0,1,by=0.01)),
         validity=function(object){
           ## Check that the proportion of defectives make sense
           if (any(is.na(object@pd)))
             return("Missing values in 'pd' not allowed")
           if (any(object@pd < 0.) | any(object@pd > 1.) )
             return("Proportion defectives must be in the range [0,1]")
         })

setClass("OChypergeom",
         representation("OC2c",
                        N="numeric",
                        pd="numeric"),
         contains="OC2c",
         prototype=list("OC2c", type="hypergeom", N=100, pd=(0:100)/100),
         validity=function(object){
           ## Check that the population size of of length 1
           if (length(object@N) > 1)
             return("Length of population size 'N' != 1")
           if (is.na(object@N))
             return("Missing value in 'N' not allowed")
           ## Check that the population size is not less than 1
           if (object@N < 1.)
             return("Population size 'N' must be at least 1")
           ## Check that the population size is non-negative
           if (object@N < sum(object@n))
             return("Total sample size must be less than population size 'N'")
           ## Check that the proportion of defectives make sense
           if (any(is.na(object@pd)))
             return("Missing value in 'pd' not allowed")
           if (any(object@pd < 0.) | any(object@pd > 1) )
             return("Proportion defectives 'pd' must be in the range [0,1]")
         })

setClass("OCpoisson",
         representation("OC2c",
                        pd="numeric"),
         contains="OC2c",
         prototype=list("OC2c", type="poisson",pd=seq(0,1,0.01)),
         validity=function(object){
           ## Check that the proportion of defectives make sense
           if (any(is.na(object@pd)))
             return("Missing values in 'pd' not allowed")
           if (any(object@pd < 0.))
             return("Rate of defects 'pd' must be non-negative")
         })

## ----------------------------------------------------------------------
## Methods to create new object and calculate P(accept)
## Only OC2c to be exported
## other functions are helpers only
## ----------------------------------------------------------------------

OC2c <- function(n,c,r=if (length(c)==1) c+1 else NULL,
               type=c("binomial","hypergeom", "poisson"), ...){
  ## Decide on what 'type' to use
  type <- match.arg(type)
  OCtype <- paste("OC",type,sep="")

  ## Create a new object of that type
  obj <- new(OCtype, n=n, c=c, r=r, type=type, ...)
  
  ## Evaluate the probability of acceptance for this type and given
  ## pd.
  ## First get the generic calculation function
##   OCtype <- get(paste("calc.",OCtype,sep=""))
  OCtype <- getFromNamespace(paste("calc.",OCtype,sep=""),
                ns="AcceptanceSampling")

  ## now, based on the type, decide on what to pass to the function
  ## Only need to check for existing type since new() would have stuffed up
  ## if we don't have a class for the type.
  if (type =="binomial")
    obj@paccept <- OCtype(n=obj@n, c=obj@c, r=obj@r, pd=obj@pd) 
  if (type =="hypergeom")
    obj@paccept <- OCtype(n=obj@n, c=obj@c, r=obj@r, N=obj@N, D=obj@pd*obj@N) 
  if (type =="poisson")
    obj@paccept <- OCtype(n=obj@n, c=obj@c, r=obj@r, pd=obj@pd) 

  obj
}




calc.OCbinomial <- function(n,c,r,pd)
{
  p.acc <- sapply(pd, FUN=calc.OCbinomial.pdi, n=n, c=c, r=r)
  p.acc
}

calc.OCbinomial.pdi <- function(pd,n,c,r)
{
  ## This is really a helper function - it does all the work for each
  ## value of pd.
  k.s <- length(n) ## number of stages in this sampling

  prob.acc <- function(x, n, p){
    k <- length(x)
    k1 <- k-1
    prod(dbinom(x[1:k1], n[1:k1], p))*pbinom(x[k], n[k], p)
  }

  
  for (k in 1:k.s) {
    ## For each stage, find out all the possibilities which could
    ## lead to still not having made a decision and then calculate
    ## the appropriate probabilities.

    if(k==1) {
      ## Only a single sampling stage to do - this is simple
      p.acc <- sapply(pd, FUN=function(el){
        pbinom(q=c[1],size=n[1],prob=el)})
      ## p.acc now exists and can be used in the following stages.
    }
    else if (k==2) {
      ## Two sampling stages. Needs to be handled separately from
      ## more stages due to matrix dimensions
      c.s <- c+1 ## Use to calculate limits
      r.s <- r-1 ## Use to calculate limits

      ## The possibilities which lead to a decision to be made at
      ## the second stage
      x <- data.frame(X1=seq(c.s[1], r.s[1], by=1),
                      X.last=c[2]-seq(c.s[1], r.s[1], by=1))
      p.acc <- p.acc + sum(apply(x, 1, FUN=prob.acc, n=n, p=pd))
    }
    else {
      ## More than two sampling stages.
      ## Things are more tricky.
      c.s <- c+1 ## Use to calculate limits
      r.s <- r-1 ## Use to calculate limits
      
      expand.call <- "expand.grid(c.s[k-1]:r.s[k-1]"
      for(i in 2:(k-1)){
        expand.call <- paste(expand.call,paste("c.s[k-",i,"]:r.s[k-",i,"]",sep=""),sep=",")
      }
      expand.call <- paste(expand.call,")",sep="")
      x <- eval(parse(text=expand.call)[[1]])
      x <- x[,(k-1):1]
      names(x) <- paste("X",1:(k-1),sep="")

      for(i in ncol(x):2){
        x[,i] <- x[,i]-x[,i-1]
      }
      x <- cbind(x, X.last=c[k] - rowSums(x[,1:(k-1)]))
      p.acc <- p.acc + sum(apply(x, 1, FUN=prob.acc, n=n, p=pd))
    }
  }
  return(p.acc)
}



calc.OChypergeom <- function(n,c,r,N,D)
{
  p.acc <- sapply(D, FUN=calc.OChypergeom.pdi, n=n, c=c, r=r, N=N)
  p.acc
}


## phyper(q=0,m=5,n=100-5,k=13) +
##   dhyper(x=1,m=5,n=100-5,k=13)*phyper(q=0,m=4,n=100-13-4,k=13)

calc.OChypergeom.pdi <- function(D,n,c,r,N)
{
  ## This is really a helper function - it does all the work for each
  ## value of pd.
  k.s <- length(n) ## number of stages in this sampling
  
  prob.acc <- function(x, n, N, D){
    k <- length(x) ## Number of sampling stages
    k1 <- k-1
    ## Total number of defects and total sample size taken so far.
    ## Note that 0 is prepended to indicate that at stage 1, zero
    ## defects have been found.
    x.cum <- cumsum(x)
    n.cum <- cumsum(n)
    N.cum <- N-c(0,n.cum[1:k1])
    D.cum <- D-c(0,x.cum[1:k1])

    prod(dhyper(x=x[1:k1], m=pmax(D.cum[1:k1],0),
                n=N.cum[1:k1]-pmax(D.cum[1:k1],0), k=n[1:k1]))*
      phyper(q=x[k], m=pmax(D.cum[k],0), n=N.cum[k]-pmax(D.cum[k],0), k=n[k])
  }

  
  for (k in 1:k.s) {
    ## For each stage, find out all the possibilities which could
    ## lead to still not having made a decision and then calculate
    ## the appropriate probabilities.

    if(k==1) {
      ## Only a single sampling stage to do - this is simple
      p.acc <- sapply(D, FUN=function(el){
        phyper(q=c[1], m=el, n=N-el, k=n[1])})
      ## p.acc now exists and can be used in the following stages.
    }
    else if (k==2) {
      ## Two sampling stages. Needs to be handled separately from
      ## more stages due to matrix dimensions
      c.s <- c+1 ## Use to calculate limits
      r.s <- r-1 ## Use to calculate limits

      ## The possibilities which lead to a decision to be made at
      ## the second stage
      x <- data.frame(X1=seq(c.s[1], r.s[1], by=1),
                      X.last=c[2]-seq(c.s[1], r.s[1], by=1))
      p.acc <- p.acc + sum(apply(x, 1, FUN=prob.acc, n=n, N=N, D=D))
    }
    else {
      ## More than two sampling stages.
      ## Things are more tricky.
      c.s <- c+1 ## Use to calculate limits
      r.s <- r-1 ## Use to calculate limits
      
      expand.call <- "expand.grid(c.s[k-1]:r.s[k-1]"
      for(i in 2:(k-1)){
        expand.call <- paste(expand.call,paste("c.s[k-",i,"]:r.s[k-",i,"]",sep=""),sep=",")
      }
      expand.call <- paste(expand.call,")",sep="")
      x <- eval(parse(text=expand.call)[[1]])
      x <- x[,(k-1):1]
      names(x) <- paste("X",1:(k-1),sep="")

      for(i in ncol(x):2){
        x[,i] <- x[,i]-x[,i-1]
      }
      x <- cbind(x, X.last=c[k] - rowSums(x[,1:(k-1)]))
      p.acc <- p.acc + sum(apply(x, 1, FUN=prob.acc, n=n, N=N, D=D))
    }
  }
  return(p.acc)
}


calc.OCpoisson <- function(n,c,r,pd)
{
  p.acc <- sapply(pd, FUN=calc.OCpoisson.pdi, n=n, c=c, r=r)
  p.acc
}

## ppois(q=c, lambda=el*n) el=pd.

calc.OCpoisson.pdi <- function(pd,n,c,r)
{
  ## This is really a helper function - it does all the work for each
  ## value of pd.
  k.s <- length(n) ## number of stages in this sampling

  prob.acc <- function(x, n, p){
    k <- length(x)
    k1 <- k-1
    prod(dpois(x[1:k1], n[1:k1]*p))*ppois(x[k], n[k]*p)
  }

  
  for (k in 1:k.s) {
    ## For each stage, find out all the possibilities which could
    ## lead to still not having made a decision and then calculate
    ## the appropriate probabilities.

    if(k==1) {
      ## Only a single sampling stage to do - this is simple
      p.acc <- sapply(pd, FUN=function(el){
        ppois(q=c[1],lambda=n[1]*el)})
      ## p.acc now exists and can be used in the following stages.
    }
    else if (k==2) {
      ## Two sampling stages. Needs to be handled separately from
      ## more stages due to matrix dimensions
      c.s <- c+1 ## Use to calculate limits
      r.s <- r-1 ## Use to calculate limits

      ## The possibilities which lead to a decision to be made at
      ## the second stage
      x <- data.frame(X1=seq(c.s[1], r.s[1], by=1),
                      X.last=c[2]-seq(c.s[1], r.s[1], by=1))
      p.acc <- p.acc + sum(apply(x, 1, FUN=prob.acc, n=n, p=pd))
    }
    else {
      ## More than two sampling stages.
      ## Things are more tricky.
      c.s <- c+1 ## Use to calculate limits
      r.s <- r-1 ## Use to calculate limits
      
      expand.call <- "expand.grid(c.s[k-1]:r.s[k-1]"
      for(i in 2:(k-1)){
        expand.call <- paste(expand.call,paste("c.s[k-",i,"]:r.s[k-",i,"]",sep=""),sep=",")
      }
      expand.call <- paste(expand.call,")",sep="")
      x <- eval(parse(text=expand.call)[[1]])
      x <- x[,(k-1):1]
      names(x) <- paste("X",1:(k-1),sep="")

      for(i in ncol(x):2){
        x[,i] <- x[,i]-x[,i-1]
      }
      x <- cbind(x, X.last=c[k] - rowSums(x[,1:(k-1)]))
      p.acc <- p.acc + sum(apply(x, 1, FUN=prob.acc, n=n, p=pd))
    }
  }
  return(p.acc)
}


## calc.OCpoisson <- function(n,c,r,pd)
## {
##   ## n needs to be cumulative since c and r are specified that way too.
##   n <- cumsum(n)
##   ## Get a list with a vector for each pd.
##   ## Length of vector equals number of samples, e.g. double = length 2.
##   ## The rate of defects is given per item.  Need to convert to
##   ## rate per sample size (multiply by n)
##   p.accept <- lapply(pd, FUN=function(el) ppois(q=c, lambda=el*n) )
##   p.unsure <- lapply(pd, FUN=function(el) {
##     ppois(q=(r-1), lambda=el*n) - ppois(q=c, lambda=el*n)})

##   ## Now combine the sampling stages via helper function
##   pa <- mapply(FUN=calc.paccept, p.accept=p.accept, p.unsure=p.unsure)
##   pa
## }

        

## ----------------------------------------------------------------------
## Printing methods and functions
## ----------------------------------------------------------------------

OC2c.show.default <-
  function(object){
    if(length(object@n)==0){
      x <- matrix(rep(NA,3), ncol=1)
    }
    else
      x <- rbind(object@n, object@c, object@r)
    dimnames(x) <- list(c("Sample size(s)", "Acc. Number(s)",
                          "Rej. Number(s)"),
                        paste("Sample", 1:ncol(x)))
    show(x)
  }

OC2c.show.prob <-
  function(object) {
    if (object@type=="binomial") {
      x <- cbind(object@pd, object@paccept)
      colnames(x) <- c("Prop. defective","P(accept)")
    }
    else if (object@type=="hypergeom"){
      x <- cbind(object@pd*object@N, object@pd, object@paccept)
      colnames(x) <- c("Pop. Defectives", "Pop. Prop. defective","P(accept)")
    }
    else if (object@type=="poisson"){
      x <- cbind(object@pd, object@paccept)
      colnames(x) <- c("Rate of defects","P(accept)")
    }
    else
      stop("No full print method defined for this type")
    
    rownames(x) <- rep("", length(object@paccept))
    show(x)
  }


setMethod("show", "OC2c",
          function(object){
            cat(paste("Acceptance Sampling Plan (",object@type,")\n\n",sep=""))
            OC2c.show.default(object)
          })

setMethod("show", "OChypergeom",
          function(object){
            cat(paste("Acceptance Sampling Plan (",
                      object@type," with N=",object@N,")\n\n",sep=""))
            OC2c.show.default(object)
          })

setMethod("summary", "OC2c",
          function(object, full=FALSE){
            cat(paste("Acceptance Sampling Plan (",object@type,")\n\n",sep=""))
            OC2c.show.default(object)
            if (full){
              cat("\nDetailed acceptance probabilities:\n\n")
              OC2c.show.prob(object)
            }
          })

setMethod("summary", "OChypergeom",
          function(object, full=FALSE){
            cat(paste("Acceptance Sampling Plan (",
                      object@type," with N=",object@N,")\n\n",sep=""))
            OC2c.show.default(object)
            if (full){
              cat("\nDetailed acceptance probabilities:\n\n")
              OC2c.show.prob(object)
            }
          })



## ----------------------------------------------------------------------
## Plotting methods
## ----------------------------------------------------------------------

setMethod("plot", signature(x="OCbinomial", y="missing"),
          function(x, y, type="o", ylim=c(0,1),...){
            plot(x@pd, x@paccept, type=type,
                 xlab="Proportion defective", ylab="P(accept)",
                 ylim=ylim, ...)
          })

setMethod("plot", signature(x="numeric", y="OCbinomial"),
          function(x, y, type="o", ylim=c(0,1),...){
            plot(x, y@paccept, type=type,
                 ylab="P(accept)", ylim=ylim, ...)
          })


setMethod("plot", signature(x="OChypergeom", y="missing"),
          function(x, type="p", ylim=c(0,1), axis=c("pd","D","both"), ...){
            xs <- match.arg(axis)

            if (xs=="pd")
              plot(x@pd, x@paccept, type=type,
                   xlab=paste("Proportion of population defectives (N=",x@N,")",sep=""),
                   ylab="P(accept)", ylim=ylim, ...)
            else if (xs=="D")
              plot(x@pd*x@N, x@paccept, type=type,
                   xlab=paste("Population defectives, D (N=",x@N,")",sep=""),
                   ylab="P(accept)", ylim=ylim, ...)
            else if (xs=="both") {
              plot(x@pd, x@paccept, type=type,
                   xlab=paste("Proportion of population defectives, (N=",x@N,")",sep=""),
                   ylab="P(accept)", ylim=ylim, mar=c(5,4,5,2)+0.1,...)
              ax <- axis(1)
              axis(3, at=ax, labels=ax*x@N)
              mtext(paste("Population defectives, D (N=",x@N,")",sep=""),
                    side=3, line=3)
            }
          })

setMethod("plot", signature(x="numeric", y="OChypergeom"),
          function(x, y, type="p", ylim=c(0,1), ...){
            plot(x, y@paccept, type=type,
                 ylab="P(accept)", ylim=ylim, ...)
          })


setMethod("plot", signature(x="OCpoisson", y="missing"),
          function(x, y, type="o", ylim=c(0,1),...){
            plot(x@pd, x@paccept, type=type,
                 xlab="Rate of defects", ylab="P(accept)",
                 ylim=ylim, ...)
          })

setMethod("plot", signature(x="numeric", y="OCpoisson"),
          function(x, y, type="o", ylim=c(0,1),...){
            plot(x, y@paccept, type=type,
                 ylab="P(accept)", ylim=ylim, ...)
          })



## ----------------------------------------------------------------------
## Methods to evaluation risk points
## All these functions are helpers only and should not be exported
## "assess" methods are exported
## ----------------------------------------------------------------------

assess.OC2c <-
  function(object, PRP, CRP){
    ## Purpose: This is the function that does the work.
    ##          Evaluate whether a particular sampling plan can meet
    ##          specified producer and/or consumer risk points
    ## ----------------------------------------------------------------------
    ## Arguments:
    ## object: An object of class OC2c
    ## PRP   : Producer risk point in the form c(pdefect, paccept)
    ## CRP   : Consumer risk point in the form c(pdefect, paccept)
    ## print : Print the result
    ## ----------------------------------------------------------------------
    ## Author: Andreas Kiermeier, Date: 16 May 2007, 10:19
    
    planOK <- TRUE
    ## Check that what we are given is OK
    if (missing(PRP))
      PRP <- rep(NA,3)
    else if (!missing(PRP)){
      if( !check.quality(PRP[1], type=object@type) |
         !check.paccept(PRP[2]) )
        stop("Quality and/or desired P(accept) out of bounds")

      ## Get the appropriate function for the distribution and
      ## calculate the P(accept)
##       calc.pa <- get(paste("calc.OC",object@type,sep=""))
      calc.pa <- getFromNamespace(paste("calc.OC",object@type,sep=""),
                     ns="AcceptanceSampling")

      pa <- switch(object@type,
                   binomial=calc.pa(object@n, object@c, object@r, PRP[1]),
                   hypergeom=calc.pa(object@n, object@c, object@r, object@N, PRP[1]*object@N),
                   poisson=calc.pa(object@n, object@c, object@r, PRP[1]))
      
      PRP <- c(PRP, pa)

      ## Check that the plan meets the desired point
      ## For PRP have to have P(accept) greater than desired prob.
      if (pa >= PRP[2])
        planOK <- TRUE
      else
        planOK <- FALSE
    }

    
    if (missing(CRP))
      CRP <- rep(NA,3)
    else if (!missing(CRP)){
      if( !check.quality(CRP[1], type=object@type) |
         !check.paccept(CRP[2]) )
        stop("Quality and/or desired P(accept) out of bound")
      ## Get the appropriate function for the distribution and
      ## calculate the P(accept)
##       calc.pa <- get(paste("calc.OC",object@type,sep=""))
      calc.pa <- getFromNamespace(paste("calc.OC",object@type,sep=""),
                     ns="AcceptanceSampling")
      pa <- switch(object@type,
                   binomial=calc.pa(object@n, object@c, object@r, CRP[1]),
                   hypergeom=calc.pa(object@n, object@c, object@r, object@N, CRP[1]*object@N),
                   poisson=calc.pa(object@n, object@c, object@r, CRP[1]))

      CRP <- c(CRP, pa)
      ## Check that the plan meets the desired point
      ## For CRP have to have P(accept) less than desired prob.
      if (pa <= CRP[2])
        planOK <- planOK & TRUE
      else
        planOK <- planOK & FALSE
    }
    return(list(OK=planOK, PRP=PRP, CRP=CRP))
  }

setMethod("assess", signature(object="OC2c"),
          function(object, PRP, CRP, print)
          {
            ## Purpose: Evaluate whether a particular sampling plan can meet
            ##          specified producer and/or consumer risk points
            ## ----------------------------------------------------------------------
            ## Arguments:
            ## object: An object of class OC2c
            ## PRP   : Producer risk point in the form c(pdefect, paccept)
            ## CRP   : Consumer risk point in the form c(pdefect, paccept)
            ## print : Print the result
            ## ----------------------------------------------------------------------
            ## Author: Andreas Kiermeier, Date: 16 May 2007, 10:19

            if(!hasArg(PRP) & !hasArg(CRP))
              stop("At least one risk point, PRP or CRP, must be specified")
            else if(CRP[1] <= PRP[1])
              stop("Consumer Risk Point quality must be greater than Producer Risk Point quality")

            plan <- assess.OC2c(object, PRP, CRP)
            if (print) {
              cat(paste("Acceptance Sampling Plan (",object@type,")\n\n",sep=""))
              OC2c.show.default(object)
              cat(paste("\nPlan", ifelse(plan$OK, "CAN","CANNOT"),
                        "meet desired risk point(s):\n\n"))

              ## Both PRP and CRP
              if(hasArg(PRP) & hasArg(CRP))
                RP <- cbind(PRP=plan$PRP, CRP=plan$CRP)
              ## Only PRP
              else if (hasArg(PRP))
                RP <- cbind(PRP=plan$PRP)
              ## Only CRP
              else if (hasArg(CRP))
                RP <- cbind(CRP=plan$CRP)

              rownames(RP) <- c("       Quality", "  RP P(accept)", "Plan P(accept)")
              show(t(RP))
            }

            if(object@type=="hypergeom")
              return(invisible(c(list(n=object@n, c=object@c, r=object@r,
                                      n=object@N), plan)))
            else
              return(invisible(c(list(n=object@n, c=object@c, r=object@r), plan)))
          })


### Local Variables:
### comment-start: "## "
### fill-column: 80
### End:
