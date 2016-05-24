# Mexico Toolkit
# Author(s) : H. Monod, INRA-MIA-Jouy, 78352, Jouy en Josas, France 
# Repository: https://mulcyber.toulouse.inra.fr/projects/baomexico/
# This file was generated automatically by the tool: mtk.AnalyserAddons() built by Juhui WANG, INRA-MIA, JOUY, FRANCE.


#' A sub-class of the class \code{\linkS4class{mtkAnalyser}} used to perform the sensitivity  analysis
#' with the "Morris" method defined in the "morrisAnalysor.R" file.
#' For more details, see the help of the function "analysor.morris()" defined in the "morrisAnalysor.R" file.
#' @title The mtkMorrisAnalyser class
#' @exportClass mtkMorrisAnalyser

setClass("mtkMorrisAnalyser",
		contains=c("mtkAnalyser")
	)

#' The constructor.
#' @param mtkParameters a vector of [\code{\linkS4class{mtkParameter}}] representing the parameters necessary to run the Analyser.
#' @param listParameters a named list defining the parameters necessary to run the Analyser. It gives non object-oriented way to specify the parameters.
#' @return an object of class \code{\linkS4class{mtkMorrisAnalyser}}
#' @examples mtkMorrisAnalyser()
#' @export mtkMorrisAnalyser
#' @title The constructor

mtkMorrisAnalyser<- function(mtkParameters=NULL, listParameters=NULL) {
	p<-mtkParameters
		if(!is.null(listParameters))
			p <- make.mtkParameterList(listParameters)
					
	res <- new("mtkMorrisAnalyser",service="Morris", parameters=p)
		return(res)
		}


#' Performs  sensitivity analysis  with the method  "Morris" defined in the "morrisAnalysor.R" file. 
#' @title The run method
#' @param this an object of class \code{\linkS4class{mtkMorrisAnalyser}}
#' @param context an object of class \code{\linkS4class{mtkExpWorkflow}}
#' @return invisible()
#' @exportMethod run
	
setMethod(f="run", signature=c(this="mtkMorrisAnalyser",
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
	
	
					analysisOutput<-eval(do.call("analysor.morris",c(list(X=X,Y=Y), parameters)))

			##!!
			##!!  post-processing the output of the method:
			##!!
	
		this@result <- mtkMorrisAnalyserResult(main=analysisOutput$main, information=analysisOutput$information)
		this@state<-TRUE
			
			assign(nameThis, this, envir=parent.frame())
			return(invisible())
		})
	
				####################### 
				## THE ORIGINAL CODE ## 
				####################### 


#' An analysor by the Morris method, using its implementation in the sensitivity library
#' @param X a design dataframe
#' @param Y a response dataframe
#' @param parameters a vector of parameters needed for the call to tell in the R sensitivity library (see details)
#' @value a list whose first element "main" is the output from the morris function and whose second element "information" is a list 
#' to give complementary information about the analysor.
#' @note The 'parameters' argument can be given the 'information' slot of the output of the 'sampler.morris' function

analysor.morris <- function(X, Y, r.design, design, scale, binf, bsup, ...){
  ##!!
  ##!! Pre-processing the input data, the processing of the method to implement follows:
  ##!!
  ## !!!! ATTENTION: AJOUTER binf et bsup !!!!
  Obj <- list(X=X, model=NULL, factors=names(X), r=r.design, design=design, scale=scale, binf=binf, bsup= bsup)
  class(Obj) <- "morris"
  
  ## Management of the multivariate case (VERY BASIC at the moment)
  Y <- Y[,1]
  
  ## Main stage:
  ## here, just a call to the "tell" function of "sensitivity"
  resultat <- list(main=tell(Obj, Y), information=NULL)
  
  ## Output
  return(resultat)
	
}

summary.mtknewmorris <- function(Object,...){
  print(Object$main)
}

plot.mtknewmorris <- function(x,y, ...){
  plot(x$main)
}
# Mexico Toolkit
# Author(s) : H. Monod, INRA-MIA-Jouy, 78352, Jouy en Josas, France 
# Repository: https://mulcyber.toulouse.inra.fr/projects/baomexico/
# This file was generated automatically by the tool: mtk.AnalyserAddons() built by Juhui WANG, INRA-MIA, JOUY, France.


#' A sub-class of the class \code{\linkS4class{mtkAnalyserResult}} used to hold the results of the sensitivity  analysis
#' with the "Morris" method defined in the "morrisAnalysor.R" file.
#' For more details, see the help of the function "analysor.morris()" defined in the "morrisAnalysor.R" file.
#' @title The mtkMorrisAnalyserResult class
#' @exportClass mtkMorrisAnalyserResult

setClass("mtkMorrisAnalyserResult",
	
				contains=c("mtkAnalyserResult")
		)

#' The constructor.
#'  @param main a data-frame to hold the main results produced by the Analyser.
#'  @param information a named list to provide supplementary information about the analysis process and its results.
				
#' @return an object of class \code{\linkS4class{mtkMorrisAnalyserResult}}
#' @examples mtkmtkMorrisAnalyserResult()
#' @export mtkmtkMorrisAnalyserResult
#' @title The constructor

mtkMorrisAnalyserResult <- function(main, information=NULL) {
	res <- new("mtkMorrisAnalyserResult", main=main, information=information)
				return(res)
				}
#' Shows a summary of the results  produced by the Analyser.
#' For more details, see the help of the function "summary.mtknewmorris()" defined in the "morrisAnalysor.R" file.
#' @title The summary method
#' @param object an object of class \code{\linkS4class{mtkMorrisAnalyserResult}}
#' @return invisible()
#' @exportMethod summary

setMethod(f="summary", "mtkMorrisAnalyserResult",
	definition=function(object,...){
summary.mtknewmorris(list(main=object@main,information=object@information),...)
	})

#' Plots the results  produced by the Analyser.
#' For more details, see the help of the function "plot.mtknewmorris()" defined in the "morrisAnalysor.R" file.
#' @title The plot method
#' @param x an object of class \code{\linkS4class{mtkMorrisAnalyserResult}}

#' @return invisible()
#' @exportMethod plot

setMethod(f="plot", "mtkMorrisAnalyserResult",
	definition=function(x,y,...){
if(!missing(y)) plot.mtknewmorris(list(main=x@main,information=x@information),y, ...)
else plot.mtknewmorris(list(main=x@main,information=x@information), ...)
	})
