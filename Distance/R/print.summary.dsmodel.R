#' Print summary of distance detection function model object
#' 
#' Provides a brief summary of a distance sampling analysis. Including: detection 
#' function parameters, model selection criterion, and optionally abundance in the
#' covered (sampled) region and its standard error.
#' 
#' @aliases print.summary.dsmodel
#' @param x a summary of distance sampling analysis
#' @param \dots unspecified and unused arguments for S3 consistency
#' @return Nothing, just prints the summary.
#' @author David L. Miller and Jeff Laake
#' @seealso \code{\link{summary.ds}}
#' @keywords utility
#' @export
print.summary.dsmodel <- function (x,...){

  # split up the object
  dht.obj <- x$dht
  model <- x$ddf
  x <- x$ds

  # routine from dht to print the tables...
  # stolen from print.dht, removed vcmatrices and cor arguments
  print.tables <- function(x,bysample){
     cat("\nSummary statistics:\n")
     print(x$summary)
     if("N" %in% names(x))
     {
        cat("\nAbundance:\n")
        print(x$N)
     }
     cat("\nDensity:\n")
     print(x$D)
     if(bysample)
     {
        cat("\nEstimates by sample:\n")
        print(x$bysample)
     }
  }

  cat("\nSummary for distance analysis \n")
  cat("Number of observations : ", x$n,"\n")
  cat("Distance range         : ", x$left, " - ",x$width,"\n")

  
  cat("\nModel :",model.description(model),"\n")
  # Remind the user that monotonicity constraints were enforced
  if(x$mono & x$mono.strict){
    cat("\nStrict monotonicity constraints were enforced.\n")
  }else if(x$mono){
    cat("\nMonotonicity constraints were enforced.\n")
  }
  cat("AIC   :", x$aic, "\n")

  # parameter summaries
  cat("\nDetection function parameters\n")
  cat("Scale Coefficients: ", "\n")
  print(x$coeff$key.scale)
  if(x$key %in% c("gamma","hr")) {
    cat("\nShape parameters: ", "\n")
    print(x$coeff$key.shape)
  }
  if(!is.null(x$coeff$adj.parm)) {
    cat("\nAdjustment term parameter(s): ", "\n")
    print(x$coeff$adj.parm)
  }
  cat("\n")
  if(!is.null(x$Nhat)){
    parameters=data.frame(Estimate=c(x$average.p,x$Nhat))
    row.names(parameters)=c("Average p", "N in covered region")
    if(!is.null(x$average.p.se)){
      parameters$SE=c(x$average.p.se,x$Nhat.se)
      parameters$CV=parameters$SE/parameters$Estimate
    }
  }else{
    parameters=data.frame(Estimate=c(x$average.p))
    row.names(parameters)=c("Average p")
    if(!is.null(x$average.p.se)){
      parameters$SE=c(x$average.p.se)
      parameters$CV=parameters$SE/parameters$Estimate
    }
  }
  print(parameters)


  ## dht stuff
  x<-dht.obj
  if(!is.null(x)){
    bysample<-FALSE # for now!
    if(is.null(x$clusters)){
      print.tables(x$individuals,bysample)
    }else{
      cat("\nSummary for clusters\n")
      print.tables(x$clusters,bysample)
      cat("\nSummary for individuals\n")
      print.tables(x$individuals,bysample)
      cat("\nExpected cluster size\n")
      # Added CV as an output LJT 14/10/09
      S<-x$Expected.S
      if(!is.null(S$se.Expected.S)){
        S$cv.Expected.S<-S$se.Expected.S/S$Expected.S
        S$cv.Expected.S[S$Expected.S==0]<-0
      }
    
      print(S)
    }
  }

  invisible()
}
