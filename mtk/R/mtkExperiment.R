# Mexico Toolkit
#
# version 	: 0.01
# date		: 30 nov 2009
# MAJ   	: 8 feb 2011
# licence	: GPL

# Author(s) :Juhui Wang, MIA-Jouy en Josas, INRA, 78352
# repository: https://mulcyber.toulouse.inra.fr/projects/baomexico/
# Web		: http://www.reseau-mexico.fr/
#
#' $Rev:: 219                 $: revision number of the last spread
#' $Author:: jwang            $: author of the last spread
#' $Date:: 2011-11-03 16:11:4#$: date of the last spread

#-----------------------------------------------------------------------

##########################################################################
# mtkExperiment class,  definition and methods
#
###########################################################################
DEBUG.mtkExperiment=FALSE



#' The class mtkExperiment is a sub-class of the class \code{\linkS4class{mtkExpWorkflow}}. 
#' It provides more facilities for interactive manipulation of the workflow.
#' 
#' @title The mtkExperiment class
#' @exportClass mtkExperiment

setClass(Class="mtkExperiment",
		contains=c("mtkExpWorkflow")
)




###########################################################################
## Methods definition
###########################################################################


#' The class mtkExperiment is a  mtkWorkflow like class but more flexible to use. Different behaviors may be expected
#' by appropriately combining the parameters. For example, 1) if sampling data 
#' is produced off-line, they will be provided with the help of the parameter "XY$X" ; 2) 
#' if simulation data is produced off-line, they can be provided through the parameter "XY$Y"; 
#' 3)if the simulation data will be produced on-line with a R function, this function can be defined through the parameter "model".
#' 4) if the simulation data will be produced on-line with a class (for example,
#' named "mtkWWDMEvaluator"),  implemented within the package mtk, the parameter "model" should  contain the short name of the class (model="WWDM" for the case of the  class "mtkWWDMEvaluator" ).
#'  5) if the sampling data is produced off-line and the
#' sensitivity analysis is fulfilled with a user-defined function, the sampling data can be provided with the parameter "XY$X" and the analysis method can be 
#' specified with the help of the parameter "method" (for example, method=fast99()).
#' 6) if the sampling data is produced off-line and the
#' sensitivity analysis is fulfilled by a method implemented within the workflow, the sampling data can be provided with the parameter "XY$X" and the analysis method can be 
#' specified by assigning a string to the parameter "method" (for example, method="Fast").
#' 7)if both the sampling  and the sensitivity analysis should be fulfilled
#'  by a user defined function, this function can be specified with the parameter "method" (for example, method=fast99()). 
#' 8) if both the sampling  and the sensitivity analysis will be fulfilled by a method implemented within the "mtk" package, the method can be specified
#' by expliciting the name of the method such as method="Fast".
#' 
#' Other scenarios can be specified by right combining the parameters.
#' 
#' @param expFactors an object of class \code{\linkS4class{mtkExpFactors}}.
#' 
#' @param model  a parameter to define the model used to produce the simulation data.
#' Its value might be  "NULL", a function or a string. "NULL" deals with the case where 
#' off-line simulation data is provided with the parameter XY$Y; a fonction deals with the case where 
#' the simulation model is defined via a function R. Such function must return a named list containing at least two elements: Y (a dataframe keeping the results 
#' produced with the simulation) and information (a named list containing supplementary information used in the simulation process such as the parameters, attributes, etc.); 
#' a string deals with the case where the simulation 
#' model is implemented within the mtk package with the help of a class. For instance, if your simulation 
#' is implemented within the mtk package with the class "mtkWWDMSimultor", you just need to give the short name of the Evaluator such as model="WWDM". 
#'
#'  @param method a parameter to define the method used to fulfill the sampling or the sensitivity  analysis.
#' Like the parameter "model", its value might be  "NULL", a function or a string. "NULL" deals with the case where 
#' off-line sampling data is provided with the parameter XY$X; a fonction deals with the case where 
#' the sampling or the analysis process is defined via a function R; a string deals with the case where the sampling or 
#' the analysis process is implemented within the mtk package with the help of a class. For instance, if your Analyser
#' is implemented within the mtk package with the class "mtkFastAnalyser", you just need to give the short name of the methd such as method="Fast", etc.
#' 
#' @param XY a named list of dataframes. The element XY$X provides with the sampling data if it's produced off-line, and the 
#' element XY$Y provided with the simulation data if it's produced off-line.
#' 
#' @param information a named list containing the parametrs necessary to run the workflow. They might be concerned 
#' with the sampling, the simulation or the analysis.
#' 
#' @return an object of class \code{\linkS4class{mtkExperiment}}
#' @title The constructor
#' @export mtkExperiment

mtkExperiment= function(expFactors, design=NULL, designInfo=NULL, model=NULL, modelInfo=NULL, analyze=NULL,  analyzeInfo=NULL, XY=NULL) {
	
	
	if(is.null(analyze)) stop("analysis method missed!")
	X=NULL; Y=NULL
	if(!is.null(XY$X)) X=XY$X
	if(!is.null(XY$Y)) Y=XY$Y
	
	designer<-mtkNativeDesigner(design=design,X=X, information=designInfo)
	evaluator<-mtkNativeEvaluator(model=model,Y=Y,information=modelInfo)
	analyser<-mtkNativeAnalyser(analyze=analyze,information=analyzeInfo)
	
	res <- new("mtkExperiment",expFactors=expFactors,processesVector=c(design=designer, evaluate=evaluator,analyze=analyser))
	return(res)
}

