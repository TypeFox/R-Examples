# Mexico Toolkit
# Author(s) : H. Monod, INRA-MIA Jouy en Josas 
# Repository: https://mulcyber.toulouse.inra.fr/projects/baomexico/
# This file was generated automatically by the tool: mtk.designerAddons() built by Juhui WANG, INRA-JOUY, FRANCE.


#' A sub-class of the class \code{\linkS4class{mtkDesigner}} used to generate the experiment design
#' with the "BasicMonteCarlo" method defined in the "montecarloDesigner.R" file.
#' For more details, see the help of the function "basicMonteCarlo()" defined in the "montecarloDesigner.R" file.
#' @title The mtkBasicMonteCarloDesigner class
#' @exportClass mtkBasicMonteCarloDesigner

setClass("mtkBasicMonteCarloDesigner",
		contains=c("mtkDesigner")
	)

#' The constructor.
#' @param mtkParameters a vector of [\code{\linkS4class{mtkParameter}}] representing the parameters necessary to run the Designer.
#' @param listParameters a named list defining the parameters necessary to run the Designer. It gives non object-oriented way to define the value of the parameters.
#' @return an object of class \code{\linkS4class{mtkBasicMonteCarloDesigner}}
#' @examples mtkBasicMonteCarloDesigner()
#' @export mtkBasicMonteCarloDesigner
#' @title The constructor

mtkBasicMonteCarloDesigner <- function(mtkParameters=NULL, listParameters=NULL) {
	p<-mtkParameters
	if(!is.null(listParameters))
					p <- make.mtkParameterList(listParameters)
					
	res <- new("mtkBasicMonteCarloDesigner",service="BasicMonteCarlo", parameters=p)
	return(res)
					}


#' Generates the experiment design with the method  "BasicMonteCarlo" defined in the "montecarloDesigner.R" file
#' @title The run method
#' @param this an object of class \code{\linkS4class{mtkBasicMonteCarloDesigner}}
#' @param context an object of class \code{\linkS4class{mtkExpWorkflow}}
#' @return invisible()
#' @exportMethod run
					
setMethod(f="run", signature=c(this="mtkBasicMonteCarloDesigner",
			context="mtkExpWorkflow"),
					definition=function(this, context){
					if(this@state) return(invisible())
					nameThis<-deparse(substitute(this))
					
					# A changer pour le Designer
					## sParametres<- as.list(formals(main))[-1]
					
					expFacteurs<-context@expFactors
					parameters<-getParameters(this)
					factorNames <-getNames(expFacteurs)
					distribNames<-getDistributionNames(expFacteurs)
					distribParameters<-getDistributionParameters(expFacteurs)
					##!!
					##!! Pre-processing the input data, the processing of the main function to implement follows:
					##!!
		
	
					arg<-list(factors=factorNames,distribNames=distribNames, distribParameters=distribParameters)	
					sortie<-eval(do.call("basicMonteCarlo",c(arg, parameters)))

					##!!
					##!!  post-processing the output of the method:
					##!!
	
					this@result <- mtkBasicMonteCarloDesignerResult(main=sortie$main, information=sortie$information)
					this@state<-TRUE
					
					assign(nameThis, this, envir=parent.frame())
					return(invisible())
					})
					
				####################### 
				## THE ORIGINAL CODE ## 
				####################### 


#' A basic Monte Carlo Designer
#' @param factors a vector of factors names
#' @param distribNames the name of a probability distribution or a vector of such names
#' @param distribParameters a list of lists of distribution parameters
#' @param size the required sample size
#' @param ...
#' @value a named list consisting of a design data.frame (main) and a list of information parameters (information)
#' @note this function uses the Quantiles function defined in the mtk package (file mtkCommonFunctions.R)
#' @examples
#'   MC <- basicMonteCarlo(LETTERS[1:4])
basicMonteCarlo <- function(factors, distribNames="unif", distribParameters=list(min=0,max=1), size, ...){

  ## PRELIMINARIES
  ## number of factors and their distributions
  nbf <- length(factors)
  if(length(distribNames) == 1){
    distribNames <- rep(distribNames, nbf)
    distribParameters <- rep(list(distribParameters), nbf)
  }
  
  ## MAIN CALCULATIONS
  N <- size*nbf
  design <- matrix(runif(N), nrow=size, ncol=nbf)
  ## quantile calculations
  for(i in seq(nbf)){
    design[,i] <- Quantiles(design[,i],	distribNames[i], distribParameters[[i]])
  }
  colnames(design) <- factors
  
  ## OUTPUT
  information <- list(SamplingMethod="Basic Monte Carlo")
  resultat <- list(main=as.data.frame(design), information=information)
  
  return(resultat)
}
# Mexico Toolkit
# Author(s) : H. Monod, INRA-MIA Jouy en Josas 
# Repository: https://mulcyber.toulouse.inra.fr/projects/baomexico/
# This file was generated automatically by the tool: mtk.AnalyserAddons(), which is built by Juhui WANG, INRA-MIA, Jouy.


#' A sub-class of the class \code{\linkS4class{mtkDesignerResult}} used to hold the results of the experiment design
#' with the "BasicMonteCarlo" method defined in the "montecarloDesigner.R" file.
#' For more details, see the help of the function "basicMonteCarlo()" defined in the "montecarloDesigner.R" file.
#' @title The mtkBasicMonteCarloDesignerResult class
#' @exportClass mtkBasicMonteCarloDesignerResult

setClass("mtkBasicMonteCarloDesignerResult",
	
					contains=c("mtkDesignerResult")
		)

#' The constructor.
#'  @param main a data frame to hold the main results produced by the Designer.
#'  @param information a named list to provide supplementary information about the sampling process and its results.
					
#' @return an object of class \code{\linkS4class{mtkBasicMonteCarloDesignerResult}}
#' @examples mtkmtkBasicMonteCarloDesignerResult()
#' @export mtkmtkBasicMonteCarloDesignerResult
#' @title The constructor

mtkBasicMonteCarloDesignerResult <- function(main=NULL, information=NULL) {
	res <- new("mtkBasicMonteCarloDesignerResult", main=main, information=information)
					return(res)
					}
