## code_var.R --- 
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
## ----------------------------------------------------------------------

## --------------------------------------------------------------------------------
## Variables Sampling Plans
## --------------------------------------------------------------------------------

setClass("OCvar", representation(n="numeric", ## A vector of sample sizes at each
                                 ## stage of sampling - NOT comulative sample size
                                 k="numeric", ## vector used to determine acceptance
                                 type="character",
                                 paccept="numeric",
                              "VIRTUAL"),
         validity=function(object){
           if(any(is.na(object@n)) | any(is.na(object@k)) )
             return("Missing values in 'n' or 'k'")
           ## Check that n and k are of length 1
           if (length(object@n) != 1 | length(object@k) != 1)
             return("n and k must be of length 1.")
           ## Check that the sample size makes sense
           if (any(object@n <= 0))
             return("Sample size 'n' must be greater than 0.")
           ## Check that the value for k makes sense
           if (any(object@k <= 0))
             return("Cut-off 'k' must be greater than 0.")
           ## Otherwise things seem fine.
           return(TRUE)
         })

setClass("OCnormal",
         representation("OCvar",
                        pd="numeric",
                        s.type="character"),
         contains="OCvar",
         prototype=list("OCvar",type="normal",pd=seq(0,1,by=0.01),s.type="known"),
         validity=function(object){
##            ## Check that the standard deviation is positive
##            if (length(object@s) > 1)
##              return("Standard deviation 's' must be of length 1.")
##            ## Check that the standard deviation is positive
##            if (object@s <= 0.)
##              return("Standard deviation 's' must be greater than 0.")
           ## Check that the s.type is either 'known' or ' unknown'
           if (object@s.type != "known" & object@s.type != "unknown")
             return("s.type must be either 'known' or 'unknown'.")
           ## Check that the proportion of defectives make sense
           if (any(is.na(object@pd)))
             return("Missing values in 'pd' not allowed")
           if (any(object@pd < 0.) | any(object@pd > 1.) )
             return("Proportion defectives must be in the range [0,1]")
         })


## ----------------------------------------------------------------------
## Methods to create new object and calculate P(accept)
## other functions are helpers only
## ----------------------------------------------------------------------

OCvar <- function(n, k, type=c("normal"), ...){
  ## Decide on what 'type' to use
  type <- match.arg(type)
  OCtype <- paste("OC",type,sep="")

  ## Create a new object of that type
  obj <- new(OCtype, n=n, k=k, type=type, ...)
  
  ## Evaluate the probability of acceptance for this type and given
  ## pd.
  ## First get the generic calculation function
  ## use 'get' for development and 'getFromNamespace' for the actual
  ## package implementation
  OCtype <- getFromNamespace(paste("calc.",OCtype,sep=""),
                ns="AcceptanceSampling")
##   OCtype <- get(paste("calc.",OCtype,sep=""))

  ## now, based on the type, decide on what to pass to the function
  ## Only need to check for existing type since new() would have stuffed up
  ## if we don't have a class for the type.

  if (type =="normal")
    obj@paccept <- OCtype(n=obj@n, k=obj@k, pd=obj@pd,
                          s.type=obj@s.type) 

  obj
}

calc.OCnormal <- function(n,k,pd,s.type)
{
  ## Are we dealing with the standard normal?
  if (s.type=="known"){
    pa <- 1-pnorm( (k+qnorm(pd))*sqrt(n))
  }
  if (s.type=="unknown"){
    pa <- 1- pt(k*sqrt(n), df=n-1, ncp=-qnorm(pd)*sqrt(n)) 
  }

  return(pa)
}

        

## ----------------------------------------------------------------------
## Printing methods and functions
## ----------------------------------------------------------------------

OCvar.show.default <-
  function(object){
    if(length(object@n)==0){
      x <- matrix(rep(NA,3), ncol=1)
    }
    else
      x <- rbind(object@n, object@k)
    dimnames(x) <- list(c("Sample size",
                          "Constant k"),
                        paste("Sample", 1:ncol(x)))
    show(x)
  }

OCvar.show.prob <-
  function(object) {
    if (object@type=="normal") {
      x <- cbind(object@pd, object@paccept)
      colnames(x) <- c("Prop. defective","P(accept)")
    }
    else
      stop("No full print method defined for this type")
    
    rownames(x) <- rep("", length(object@paccept))
    show(x)
  }

setMethod("show", "OCvar",
          function(object){
            cat(paste("Acceptance Sampling Plan (",object@type,")\n",sep=""))
            cat(paste("Standard deviation assumed to be ",object@s.type,"\n\n",sep=""))
            OCvar.show.default(object)
          })

setMethod("summary", "OCvar",
          function(object, full=FALSE){
            show(object)
            if (full){
              cat("\nDetailed acceptance probabilities:\n\n")
              OCvar.show.prob(object)
            }
          })


## ----------------------------------------------------------------------
## Plotting methods
## ----------------------------------------------------------------------

setMethod("plot", signature(x="OCnormal", y="missing"),
          function(x, y, type="o", ylim=c(0,1),...){
            plot(x@pd, x@paccept, type=type,
                 xlab="Proportion defective", ylab="P(accept)",
                 ylim=ylim, ...)
          })

setMethod("plot", signature(x="numeric", y="OCnormal"),
          function(x, y, type="o", ylim=c(0,1),...){
            plot(x, y@paccept, type=type,
                 ylab="P(accept)", ylim=ylim, ...)
          })



## ----------------------------------------------------------------------
## Methods to evaluation risk points
## All these functions are helpers only and should not be exported
## "assess" methods are exported
## ----------------------------------------------------------------------

assess.OCvar <-
  function(object, PRP, CRP){
    ## Purpose: This is the function that does the work.
    ##          Evaluate whether a particular sampling plan can meet
    ##          specified producer and/or consumer risk points
    ## ----------------------------------------------------------------------
    ## Arguments:
    ## object: An object of class OCvar
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
      calc.pa <- getFromNamespace(paste("calc.OC",object@type,sep=""),
                     ns="AcceptanceSampling")
##       calc.pa <- get(paste("calc.OC",object@type,sep=""))
      
      pa <- switch(object@type,
                   normal=calc.pa(n=object@n, k=object@k, s.type=object@s.type,
                   pd=PRP[1]))
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
      calc.pa <- getFromNamespace(paste("calc.OC",object@type,sep=""),
                     ns="AcceptanceSampling")
##       calc.pa <- get(paste("calc.OC",object@type,sep=""))

      pa <- switch(object@type,
                   normal=calc.pa(n=object@n, k=object@k, s.type=object@s.type,
                     pd=CRP[1]))
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


setMethod("assess", signature(object="OCvar"),
          function(object, PRP, CRP, print)
          {
            ## Purpose: Evaluate whether a particular sampling plan can meet
            ##          specified producer and/or consumer risk points
            ## ----------------------------------------------------------------------
            ## Arguments:
            ## object: An object of class OCvar
            ## PRP   : Producer risk point in the form c(pdefect, paccept)
            ## CRP   : Consumer risk point in the form c(pdefect, paccept)
            ## print : Print the result
            ## ----------------------------------------------------------------------
            ## Author: Andreas Kiermeier, Date: 21 August 2007

            if(!hasArg(PRP) & !hasArg(CRP))
              stop("At least one risk point, PRP or CRP, must be specified")
            else if(CRP[1] <= PRP[1])
              stop("Consumer Risk Point quality must be greater than Producer Risk Point quality")

            plan <- assess.OCvar(object, PRP, CRP)
            if (print) {
              show(object)
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

            return(invisible(c(list(n=object@n, k=object@k, s.type=object@s.type), plan)))
          })







### Local Variables:
### comment-start: "## "
### fill-column: 80
### End:
