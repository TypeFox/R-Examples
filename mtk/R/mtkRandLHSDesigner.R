# Mexico Toolkit
# Author(s) : R.Faivre (2012) 
# Repository: https://mulcyber.toulouse.inra.fr/projects/baomexico/
# This file was generated automatically by the tool: mtk.designerAddons() built by Juhui WANG, INRA-JOUY, FRANCE.


#' A sub-class of the class \code{\linkS4class{mtkDesigner}} used to generate the experiment design
#' with the "RandLHS" method defined in the "RandomLHS.R" file.
#' For more details, see the help of the function "RandomLHS()" defined in the "RandomLHS.R" file.
#' @title The mtkRandLHSDesigner class
#' @exportClass mtkRandLHSDesigner

setClass("mtkRandLHSDesigner",
		contains=c("mtkDesigner")
	)

#' The constructor.
#' @param mtkParameters a vector of [\code{\linkS4class{mtkParameter}}] representing the parameters necessary to run the Designer.
#' @param listParameters a named list defining the parameters necessary to run the Designer. It gives non object-oriented way to define the value of the parameters.
#' @return an object of class \code{\linkS4class{mtkRandLHSDesigner}}
#' @examples mtkRandLHSDesigner()
#' @export mtkRandLHSDesigner
#' @title The constructor

mtkRandLHSDesigner <- function(mtkParameters=NULL, listParameters=NULL) {
	p<-mtkParameters
	if(!is.null(listParameters))
					p <- make.mtkParameterList(listParameters)
					
	res <- new("mtkRandLHSDesigner",service="RandLHS", parameters=p)
	return(res)
					}


#' Generates the experiment design with the method  "RandLHS" defined in the "RandomLHS.R" file
#' @title The run method
#' @param this an object of class \code{\linkS4class{mtkRandLHSDesigner}}
#' @param context an object of class \code{\linkS4class{mtkExpWorkflow}}
#' @return invisible()
#' @exportMethod run
					
setMethod(f="run", signature=c(this="mtkRandLHSDesigner",
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
					sortie<-eval(do.call("RandomLHS",c(arg, parameters)))

					##!!
					##!!  post-processing the output of the method:
					##!!
	
		this@result <- mtkRandLHSDesignerResult(main=sortie$main, information=sortie$information)
		this@state<-TRUE
					
					assign(nameThis, this, envir=parent.frame())
					return(invisible())
					})
					
				####################### 
				## THE ORIGINAL CODE ## 
				####################### 

RandomLHS = function(factors, distribNames, distribParameters,
size = 10, preserveDraw = FALSE, ...) 
{
## Les arguments de la fonction doivent doit être 
## factors, distribNames, distribParameters
## suivis des arguments spécifiques à la fonction encapsulée

  nbf <- length(factors)
  design = randomLHS( size, nbf, preserveDraw)
## Utilisation de distribparameters pour remettre dans les bornes

for (i in 1:nbf) {
design[,i]  = distribParameters[[i]][[1]] + design[,i] *
(distribParameters[[i]][[2]] - distribParameters[[i]][[1]])
}
  colnames(design) <- factors
## OUTPUT
  information <- list(SamplingMethod="Random LHS")
## La sortie soit être une structure liste :
## list(main = ..., information = ...)
  resultat <- list(main=as.data.frame(design), information=information)
}

# Mexico Toolkit
# Author(s) : R.Faivre (2012) 
# Repository: https://mulcyber.toulouse.inra.fr/projects/baomexico/
# This file was generated automatically by the tool: mtk.AnalyserAddons(), which is built by Juhui WANG, INRA-MIA, Jouy.


#' A sub-class of the class \code{\linkS4class{mtkDesignerResult}} used to hold the results of the experiment design
#' with the "RandLHS" method defined in the "RandomLHS.R" file.
#' For more details, see the help of the function "RandomLHS()" defined in the "RandomLHS.R" file.
#' @title The mtkRandLHSDesignerResult class
#' @exportClass mtkRandLHSDesignerResult

setClass("mtkRandLHSDesignerResult",
	
					contains=c("mtkDesignerResult")
		)

#' The constructor.
#'  @param main a data frame to hold the main results produced by the Designer.
#'  @param information a named list to provide supplementary information about the sampling process and its results.
					
#' @return an object of class \code{\linkS4class{mtkRandLHSDesignerResult}}
#' @examples mtkmtkRandLHSDesignerResult()
#' @export mtkmtkRandLHSDesignerResult
#' @title The constructor

mtkRandLHSDesignerResult <- function(main=NULL, information=NULL) {
	res <- new("mtkRandLHSDesignerResult", main=main, information=information)
					return(res)
					}
