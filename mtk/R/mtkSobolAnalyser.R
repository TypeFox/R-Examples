# Mexico Toolkit
# Author(s) : H. Monod, INRA-MIA-Jouy, 78352, Jouy en Josas, France 
# Repository: https://mulcyber.toulouse.inra.fr/projects/baomexico/
# This file was generated automatically by the tool: mtk.AnalyserAddons() built by Juhui WANG, INRA-MIA, JOUY, FRANCE.


#' A sub-class of the class \code{\linkS4class{mtkAnalyser}} used to perform the sensitivity  analysis
#' with the "Sobol" method defined in the "sobolAnalyser.R" file.
#' For more details, see the help of the function "Analyser.morris()" defined in the "sobolAnalyser.R" file.
#' @title The mtkSobolAnalyser class
#' @exportClass mtkSobolAnalyser

setClass("mtkSobolAnalyser",
		contains=c("mtkAnalyser")
	)

#' The constructor.
#' @param mtkParameters a vector of [\code{\linkS4class{mtkParameter}}] representing the parameters necessary to run the Analyser.
#' @param listParameters a named list defining the parameters necessary to run the Analyser. It gives non object-oriented way to specify the parameters.
#' @return an object of class \code{\linkS4class{mtkSobolAnalyser}}
#' @examples mtkSobolAnalyser()
#' @export mtkSobolAnalyser
#' @title The constructor

mtkSobolAnalyser<- function(mtkParameters=NULL, listParameters=NULL) {
	p<-mtkParameters
		if(!is.null(listParameters))
			p <- make.mtkParameterList(listParameters)
					
	res <- new("mtkSobolAnalyser",service="Sobol", parameters=p)
		return(res)
		}


#' Performs  sensitivity analysis  with the method  "Sobol" defined in the "sobolAnalyser.R" file. 
#' @title The run method
#' @param this an object of class \code{\linkS4class{mtkSobolAnalyser}}
#' @param context an object of class \code{\linkS4class{mtkExpWorkflow}}
#' @return invisible()
#' @exportMethod run
	
setMethod(f="run", signature=c(this="mtkSobolAnalyser",
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
	
	
					analysisOutput<-eval(do.call("Analyser.sobol",c(list(X=X,Y=Y), parameters)))

			##!!
			##!!  post-processing the output of the method:
			##!!
	
		this@result <- mtkSobolAnalyserResult(main=analysisOutput$main, information=analysisOutput$information)
		this@state<-TRUE
			
			assign(nameThis, this, envir=parent.frame())
			return(invisible())
		})
	
				####################### 
				## THE ORIGINAL CODE ## 
				####################### 


#' An Analyser by the Sobol method, using its implementation in the sensitivity library
#' @param X a design dataframe
#' @param Y a response dataframe
#' @param X1 a dataframe of the first basic design used to construct X
#' @param X2 a dataframe of the second basic design used to construct X
#' @param parameters a vector of parameters needed for the call to the tell function of the R sensitivity library (see details)
#' @param nboot see the help on function sobol2002 in the sensitivity library
#' @param conf see the help on function sobol2002 in the sensitivity library
#' @value a list whose first element "main" is the output from the sobol function, and whose second element "information" is a list 
#' to give complementary information about the Analyser.
#' @note The 'parameters' argument can be given the 'information' slot of the output of the 'Designer.sobol' function
#' @examples
#'  X <- data.frame(X1 = runif(n, 0.5, 1.5),
#'                     X2 = runif(n, 1.5, 4.5),
#'                     X3 = runif(n, 4.5, 13.5))
#'  y <- with(X, X1 + X2 + X3)
#'  Y <- regressionSI(X,data.frame(y=y, y.by.y=y*y), z=exp(z))
#'  print.regressionSI(toto)
#'  plot.regressionSI(toto)

Analyser.sobol <- function(X, Y, X1, X2, nboot=0, conf=0.95, ...){
  
  ##!!
  ##!! Pre-processing the input data, the processing of the method to implement follows:
  ##!!
  
  Obj <- list(model=NULL, X=X, X1=X1, X2=X2, SamplingMethod="sobol2002", nboot=nboot, conf=conf)
  class(Obj) <- "sobol2002"
  
  ## Management of the multivariate case (VERY BASIC at the moment)
  y <- c(Y[,1])
  
  ## Main stage:
  ## here, just a call to the "tell" function of "sensitivity"
  sobolOutput <- tell(Obj, y)
  resultat <- list(main=sobolOutput, information=NULL)
  
  ## Output
  return(resultat)
	
}

summary.sobolAnalys <- function(Object){
  cat("\nThe results of an analysis by the Sobol method\n\n")
  
  toprint <- dim(Object$main[[c("X")]])
  cat("Design size: ", toprint[1], " rows by ", toprint[2], " columns\n")

  toprint <- Object$main[c("nboot","conf")]
  cat("Parameters: nboot = ", toprint[[1]], ", conf = ", toprint[[2]], "\n")
  
  cat("\nResults:\n")
  print(Object$main[c("V","S","T")])
}

plot.sobolAnalys <- function(x,y, ...){
  plot(x$main)
}

# Mexico Toolkit
# Author(s) : H. Monod, INRA-MIA-Jouy, 78352, Jouy en Josas, France 
# Repository: https://mulcyber.toulouse.inra.fr/projects/baomexico/
# This file was generated automatically by the tool: mtk.AnalyserAddons() built by Juhui WANG, INRA-MIA, JOUY, France.


#' A sub-class of the class \code{\linkS4class{mtkAnalyserResult}} used to hold the results of the sensitivity  analysis
#' with the "Sobol" method defined in the "sobolAnalyser.R" file.
#' For more details, see the help of the function "Analyser.morris()" defined in the "sobolAnalyser.R" file.
#' @title The mtkSobolAnalyserResult class
#' @exportClass mtkSobolAnalyserResult

setClass("mtkSobolAnalyserResult",
	
				contains=c("mtkAnalyserResult")
		)

#' The constructor.
#'  @param main a data-frame to hold the main results produced by the Analyser.
#'  @param information a named list to provide supplementary information about the analysis process and its results.
				
#' @return an object of class \code{\linkS4class{mtkSobolAnalyserResult}}
#' @examples mtkmtkSobolAnalyserResult()
#' @export mtkmtkSobolAnalyserResult
#' @title The constructor

mtkSobolAnalyserResult <- function(main=NULL, information=NULL) {
	res <- new("mtkSobolAnalyserResult", main=main, information=information)
				return(res)
				}
#' Shows a summary of the results  produced by the Analyser.
#' For more details, see the help of the function "summary.sobolAnalys()" defined in the "sobolAnalyser.R" file.
#' @title The summary method
#' @param object an object of class \code{\linkS4class{mtkSobolAnalyserResult}}
#' @return invisible()
#' @exportMethod summary

setMethod(f="summary", "mtkSobolAnalyserResult",
	definition=function(object,...){
summary.sobolAnalys(list(main=object@main,information=object@information),...)
	})

#' Plots the results  produced by the Analyser.
#' For more details, see the help of the function "plot.sobolAnalys()" defined in the "sobolAnalyser.R" file.
#' @title The plot method
#' @param x an object of class \code{\linkS4class{mtkSobolAnalyserResult}}

#' @return invisible()
#' @exportMethod plot

setMethod(f="plot", "mtkSobolAnalyserResult",
	definition=function(x,y,...){
if(!missing(y)) plot.sobolAnalys(list(main=x@main,information=x@information),y, ...)
else plot.sobolAnalys(list(main=x@main,information=x@information), ...)
	})
