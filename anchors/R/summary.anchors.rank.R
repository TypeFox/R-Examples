#######################################################################
##
## Function: summary.anchors.rank()
## Author  : Jonathan Wand <wand(at)stanford.edu>
## Created :  2006-10-02
##
## 
##
## MODIFIED:
##    2007-09-01 : JW
##    - added B
## 
##    2008-04-20 : JW
##    - was summary.anchors()
##    - option to look at B and C simultaneously
##   
#######################################################################
summary.anchors.rank <- function(object, ... ,
                                 ties = c("omit","uniform","cpolr","minentropy"),
                                 combn = TRUE, digits=3) {

  ties <- match.arg(ties, c("omit","uniform","cpolr","minentropy"), several.ok=TRUE)
  
  cat("\nANCHORS: SUMMARY OF RELATIVE RANK ANALYSIS:\n")
  
  cat("\nOverview of ",object$type,"-ranks\n",sep="")
  cat("\nNumber of cases:",object$summary$n.interval,"with interval value,",object$summary$n.scalar,"with scalar value\n")
  cat("\nMaximum possible ",object$type,"-rank value: ",object$summary$max,"\n",sep="")

  if (!is.null(object$minentropy) && "minentropy"  %in% ties) {
    cat("\nInterval on ",object$type,"-scale: Frequency and proportions ",object$type,"s to ",object$type,"e\n",sep="")
    print( round( summary.minimum.entropy( object$minentropy ),digits) )

    cat("\nNote: MinEnt is the rank for the interval that minimizes entropy\n")

  } else {
    cat("\nInterval on ",object$type,"-scale: Frequency and proportions ",object$type,"s to ",object$type,"e\n",sep="")
    print(format(object$summary$interval[, -c(3:4)] ,digits=digits))
    cat("\n")
  }
  
  cat("\nSummary of ",object$type,"-ranks with ties/intervals broken:\n",sep="")

  RA <- list()
  
  if ("omit" %in% ties) {
    cat("\nDistribution of ranks omiting interval cases\n")
    out <- as.data.frame(matrix(object$summary$scalar$Prop, nrow=1))
    colnames(out) <- object$summary$scalar$from
    rownames(out) <- ""
    print(round( out,digits))
    RA$omit <- RV <- as.matrix(out)
  }
  if ("uniform"  %in% ties) {
    cat("\nDistribution of ranks allocating interval cases uniformly\n")
    out <- as.data.frame(matrix(object$summary$uniform$Prop, nrow=1))
    colnames(out) <- object$summary$uniform$from
    rownames(out) <- ""
    print(round( out,digits))
    RA$uniform <- RV <- as.matrix(out)
  }
  if (!is.null(object$cpolr) && "cpolr"  %in% ties) {
    cat("\nDistribution of ranks allocating interval cases via cpolr\n")
    cat("and conditioning on observed ranks\n")
    out <- fitted( object$cpolr, object$rank, average=TRUE, unconditional=FALSE)
    class(out) <- NULL
    print( round( out ,digits))
#    cat("\nCPOLR model parameters:\n")
#    print(summary(object$cpolr))
    RA$cpolr <- RV <- as.matrix(out)
  }
  if (is.null(object$cpolr) && "cpolr"  %in% ties) {
    cat("\nThere is no cpolr model included in anchors.rank object\n")
    cat("Either cpolr was not requested when the ranks were calculated\n")
    cat("[ see 'anchors.options(rank.extras)' ],\n")
    cat("or there were not enough different ranks (<=2) to estimate the cpolr model\n")
  }
  if (!is.null(object$minentropy) && "minentropy"  %in% ties) {
    cat("\nAllocating cases to their MinEnt values produces\n")
    print( round( out <- summary.minimum.entropy( object$minentropy , average=TRUE), digits) )
    RA$minentropy <- RV <- as.matrix(out)
  }
  
#   if (is.null(object$minentropy) && "minentropy"  %in% ties) {
#     cat("\nThere is no minentropy calcuation included in anchors.rank object\n")
#     cat("'minentropy may not have been requested when the ranks were calculated\n")
#     cat("[ see 'anchors.options(rank.extras)'] \n")
#   }

  ## combn analysis 
  sort <- "max"
  if ( !is.logical(combn) ) {
    choices <- c("max","estimated","minimum","interval","span")
    i <- pmatch( combn,  choices)
    if (any(is.na(i)) || length(i) > 1)
      stop("if 'combn' is not a logical, then should be one of ", paste(choices, collapse = ", "))
    combn <- TRUE
    sort <- choices[i]
  }

  if ( combn && !is.null(object$combn)) {
    summary.anchors.combn( object$combn)
  }

  class(RV) <- class(RA) <- "summary.anchors.rank"
  
  if (length(ties)==1) {
    return(invisible(RV))
  } else {
    return(invisible(RA))
  }
}

