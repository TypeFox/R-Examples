# Mexico Toolkit
# Author(s) : H. Monod, INRA-MIA-Jouy, 78352, Jouy en Josas, France 
# Repository: https://mulcyber.toulouse.inra.fr/projects/baomexico/
# This file was generated automatically by the tool: mtk.AnalyserAddons() built by Juhui WANG, INRA-MIA, JOUY, FRANCE.


#' A sub-class of the class \code{\linkS4class{mtkAnalyser}} used to perform the sensitivity  analysis
#' with the "Regression" method defined in the "srcAnalyser.R" file.
#' For more details, see the help of the function "regressionSI()" defined in the "srcAnalyser.R" file.
#' @title The mtkRegressionAnalyser class
#' @exportClass mtkRegressionAnalyser

setClass("mtkRegressionAnalyser",
		contains=c("mtkAnalyser")
	)

#' The constructor.
#' @param mtkParameters a vector of [\code{\linkS4class{mtkParameter}}] representing the parameters necessary to run the Analyser.
#' @param listParameters a named list defining the parameters necessary to run the Analyser. It gives non object-oriented way to specify the parameters.
#' @return an object of class \code{\linkS4class{mtkRegressionAnalyser}}
#' @examples mtkRegressionAnalyser()
#' @export mtkRegressionAnalyser
#' @title The constructor

mtkRegressionAnalyser<- function(mtkParameters=NULL, listParameters=NULL) {
	p<-mtkParameters
		if(!is.null(listParameters))
			p <- make.mtkParameterList(listParameters)
					
	res <- new("mtkRegressionAnalyser",service="Regression", parameters=p)
		return(res)
		}


#' Performs  sensitivity analysis  with the method  "Regression" defined in the "srcAnalyser.R" file. 
#' @title The run method
#' @param this an object of class \code{\linkS4class{mtkRegressionAnalyser}}
#' @param context an object of class \code{\linkS4class{mtkExpWorkflow}}
#' @return invisible()
#' @exportMethod run
	
setMethod(f="run", signature=c(this="mtkRegressionAnalyser",
			context="mtkExpWorkflow"),
		definition=function(this, context){
			if(this@state) return(invisible())
			nameThis<-deparse(substitute(this))
			
			X <- context@processesVector$design@result@main 
			Y <- context@processesVector$evaluate@result@main
			parameters<-getParameters(this)
			parameters<-c(context@processesVector$design@result@information,parameters)
			##!!
			##!! Pre-processing the input data, the processing of the method to implement follows:
			##!!
	
	
					analysisOutput<-eval(do.call("regressionSI",c(list(X=X,Y=Y), parameters)))

			##!!
			##!!  post-processing the output of the method:
			##!!
	
		this@result <- mtkRegressionAnalyserResult(main=analysisOutput$main, information=analysisOutput$information)
		this@state<-TRUE
			
			assign(nameThis, this, envir=parent.frame())
			return(invisible())
		})
	
				####################### 
				## THE ORIGINAL CODE ## 
				####################### 


#' Sensitivity indices based on simple and standard regression analyses
#' @param X a design dataframe
#' @param Y a dataframe of responses
#' @param model the regression model
#' @value a list of regression results
#' @note This function is essentially a (possibly repeated) call to the "src" function of the R sensitivity package

regressionSI <- function(X, Y, rank=FALSE, nboot=100, conf=0.95, ...){

  ## PRELIMINARIES
  nf <- ncol(Y)
  Ynames <- names(Y)
  results <- vector(mode="list", length=nf)

  ## MAIN CALCULATIONS
  for(i in seq(nf)){
    resultsA <- src(X, Y[,i], rank=rank, nboot=nboot, conf=conf)
    if(!rank)
      results[[i]] <- resultsA[c("SRC")]
    else
      results[[i]] <- resultsA[c("SRRC")]
  }
  
  ## Output
  output <- list(main=results, information=list(Ynames=Ynames,rank=rank,nboot=nboot,conf=conf))
  return(output)
}

print.regressionSI <- function(x, ...){
  ## PRELIMINARIES
  results <- x$main
  Ynames <- x$information$Ynames
  rank <- x$information$rank
  nf <- length(results)

  ## MAIN CALCULATIONS
  cat("\nStandardized Regression Coefficients (SRC)\n")
  cat(" calculated by the src function of the sensitivity library\n")
  cat(" with parameters rank = ", x$information$rank,
      ", nboot = ", x$information$nboot,
      ", conf = ", x$information$conf,
      "\n\n")
  for(i in seq(nf)){
    cat("Response variable: ", Ynames[i], "\n") 
    print(results[[i]])
  }
  invisible()
}

plot.regressionSI <- function(x,y, ...){
  ## PRELIMINARIES
  results <- x$main
  Ynames <- x$information$Ynames
  rank <- x$information$rank
  nf <- length(results)
 
  ## MAIN CALCULATIONS
  nb.col <- ceiling(sqrt(nf))
  nb.row <- ceiling(nf/nb.col)
  keep.mfrow <- par(mfrow=c(nb.row,nb.col))

  for(i in seq(nf)){
    src.obj.i <- results[[i]]
    class(src.obj.i) <- "src"
    plot(src.obj.i, ...)
    title(xlab=Ynames[i])
    abline(h=0)
  }
  par(keep.mfrow)
  
  invisible()
}

# Mexico Toolkit
# Author(s) : H. Monod, INRA-MIA-Jouy, 78352, Jouy en Josas, France 
# Repository: https://mulcyber.toulouse.inra.fr/projects/baomexico/
# This file was generated automatically by the tool: mtk.AnalyserAddons() built by Juhui WANG, INRA-MIA, JOUY, France.


#' A sub-class of the class \code{\linkS4class{mtkAnalyserResult}} used to hold the results of the sensitivity  analysis
#' with the "Regression" method defined in the "srcAnalyser.R" file.
#' For more details, see the help of the function "regressionSI()" defined in the "srcAnalyser.R" file.
#' @title The mtkRegressionAnalyserResult class
#' @exportClass mtkRegressionAnalyserResult

setClass("mtkRegressionAnalyserResult",
	
				contains=c("mtkAnalyserResult")
		)

#' The constructor.
#'  @param main a data-frame to hold the main results produced by the Analyser.
#'  @param information a named list to provide supplementary information about the analysis process and its results.
				
#' @return an object of class \code{\linkS4class{mtkRegressionAnalyserResult}}
#' @examples mtkmtkRegressionAnalyserResult()
#' @export mtkmtkRegressionAnalyserResult
#' @title The constructor

mtkRegressionAnalyserResult <- function(main=NULL, information=NULL) {
	res <- new("mtkRegressionAnalyserResult", main=main, information=information)
				return(res)
				}
#' Plots the results  produced by the Analyser.
#' For more details, see the help of the function "plot.regressionSI()" defined in the "srcAnalyser.R" file.
#' @title The plot method
#' @param x an object of class \code{\linkS4class{mtkRegressionAnalyserResult}}

#' @return invisible()
#' @exportMethod plot

setMethod(f="plot", "mtkRegressionAnalyserResult",
	definition=function(x,y,...){
if(!missing(y)) plot.regressionSI(list(main=x@main,information=x@information),y, ...)
else plot.regressionSI(list(main=x@main,information=x@information), ...)
	})
#' Prints the results  produced by the Analyser.
#' For more details, see the help of the function "print.regressionSI()" defined in the "srcAnalyser.R" file.
#' @title The print method
#' @param object an object of class \code{\linkS4class{mtkRegressionAnalyserResult}}
#' @return invisible()
#' @exportMethod print

setMethod(f="print", "mtkRegressionAnalyserResult",
	definition=function(x, ...){
print.regressionSI(list(main=x@main,information=x@information), ...)
	})
