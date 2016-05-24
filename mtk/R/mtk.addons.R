# Mexico Toolkit
#
# version 	: 0.01
# date		: 22 sept. 2011
# MAJ   	: 19 oct. 2011
# licence	: GPL


#=======================================================================
#'
#' 
#'  A function used to extend the MTK package  with existing sampling methods programmed as R functions. The Designer must have the following syntax:
#'  main <- function(factors, distribNames, distribParameters, ...)
#'  where "factors" is a number or a list  managing the names of the factors, and "distribNames" is a list managing the names
#'  of the distributions of the factor, and "distribParameters" is a list of parameters associated with the distributions.
#' 
#'  The function returns a named list with two elements: the element "main" is a data frame holding the results of the experiment design and the element "information"
#'  is a named list holding complementary information about the design.
#' 
#' 	Furthermore, users can redefine the generic functions:  summary(object, ...), plot(x,y, ...), print(x, ...).
#'	
#'	For example, a Sobol Designer programmed by H. Monod is held in file "sobolDesigner.R" where the main function is called "Designer.sobol()". We use the following call
#'	to generate the mtk compliant S4 classes:
#'
#' 	mtk.designerAddons(where="DesignerSobol.R", authors="H. Monod, INRA-MIA Jouy en Josas", name="HMSobol", main="Designer.sobol")
#'
#'	The generated file named "mtkHMSobolAnalyser.R" can be directly integrated into the mtk package.

mtk.designerAddons = function(where=NULL, library=NULL, authors=NULL, name=NULL, main=NULL,summary=NULL,plot=NULL, print=NULL)
{
	if(is.null(library)&&!file.exists(where))stop(paste("The file: ", where," or the library: ",library, "does not exist in this context!", sep=" "))
	
	className<-paste("mtk",name,"Designer", sep="")
	outFile<-file(paste(className,".R", sep=""),open="w")
	cat(file=outFile,  "# Mexico Toolkit\n")
	cat(file=outFile, "# Author(s) :", authors,"\n")
	cat(file=outFile, "# Repository: https://mulcyber.toulouse.inra.fr/projects/baomexico/\n")
	cat(file=outFile, "# This file was generated automatically by the tool: mtk.designerAddons() built by Juhui WANG, INRA-JOUY, FRANCE.\n\n\n")
	
	if(!is.null(library))cat(file=outFile, "library(",library,")\n\n")
	
	cat(file=outFile, "#' A sub-class of the class \\code{\\linkS4class{mtkDesigner}} used to generate the experiment design\n")
	if(!is.null(where))cat(file=outFile, "#' with the \"", name, "\" method defined in the \"", where, "\" file.\n", sep="")
	if(!is.null(library)) cat(file=outFile, "#' with the \"", name, "\"  method defined in the \"", library, "\" library.\n", sep="")
	if(!is.null(where))cat(file=outFile, "#' For more details, see the help of the function \"", main, "()\" defined in the \"", where, "\" file.\n", sep="")
	if(!is.null(library)) cat(file=outFile, "#' For more details, see the help of the function \"", main, "()\" defined in the \"", library, "\" library.\n", sep="")
	
	cat(file=outFile, "#' @title The mtk", name, "Designer class\n", sep="")
	cat(file=outFile, "#' @exportClass mtk",name,"Designer\n\n", sep="")
	
	cat(file=outFile, "setClass(\"mtk",name,"Designer\",\n", sep="")
	cat(file=outFile, "		contains=c(\"mtkDesigner\")\n")
	cat(file=outFile, "	)")
	
	cat(file=outFile, "\n\n#' The constructor.
#' @param mtkParameters a vector of [\\code{\\linkS4class{mtkParameter}}] representing the parameters necessary to run the Designer.
#' @param listParameters a named list defining the parameters necessary to run the Designer. It gives non object-oriented way to define the value of the parameters.\n")
	
	cat(file=outFile, "#' @return an object of class \\code{\\linkS4class{mtk",name,"Designer}}\n",sep="")
	cat(file=outFile, "#' @examples mtk",name,"Designer()\n",sep="")
	cat(file=outFile, "#' @export mtk",name,"Designer\n",sep="")
	cat(file=outFile, "#' @title The constructor\n\n")
	
	cat(file=outFile, className," <- function(mtkParameters=NULL, listParameters=NULL) {\n",sep="")
	cat(file=outFile, "	p<-mtkParameters
	if(!is.null(listParameters))
					p <- make.mtkParameterList(listParameters)
					\n")
	cat(file=outFile, "	res <- new(\"",className,"\",service=\"",name,"\", parameters=p)
	return(res)
					}\n", sep="")
	
	if(!is.null(where))cat(file=outFile, "\n\n#' Generates the experiment design with the method  \"",name,"\" defined in the \"",where,"\" file\n", sep="")
	if(!is.null(library)) cat(file=outFile, "\n\n#' Generates the experiment design with the method  \"",name,"\" defined in the \"",library,"\" library\n", sep="")
	
	cat(file=outFile, "#' @title The run method\n")
	cat(file=outFile, "#' @param this an object of class \\code{\\linkS4class{mtk", name,"Designer}}\n",sep="")
	cat(file=outFile, "#' @param context an object of class \\code{\\linkS4class{mtkExpWorkflow}}
#' @return invisible()
#' @exportMethod run
					\n")
	
	cat(file=outFile, "setMethod(f=\"run\", signature=c(this=\"", className, "\",\n", sep="")
	cat(file=outFile, "			context=\"mtkExpWorkflow\"),
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
		\n")			
	cat(file=outFile, "	
					arg<-list(factors=factorNames,distribNames=distribNames, distribParameters=distribParameters)")
		
	cat(file=outFile, "	
					sortie<-eval(do.call(\"",main,"\",c(arg, parameters)))\n",sep="")
	cat(file=outFile, "
					##!!
					##!!  post-processing the output of the method:
					##!!
	\n")
	cat(file=outFile, "		this@result <- mtk",name,"DesignerResult(main=sortie$main, information=sortie$information)\n",sep="")
	cat(file=outFile, "		this@state<-TRUE
					
					assign(nameThis, this, envir=parent.frame())
					return(invisible())
					})
					\n")
	if(!is.null(where))
	{
		cat(file=outFile, "				####################### \n")
		cat(file=outFile, "				## THE ORIGINAL CODE ## \n")
		cat(file=outFile, "				####################### \n\n\n")
		
		inFile<-file(where, open="r")
		src<-readLines(con=inFile, n=-1)
		close(outFile)
		outFile<-file(paste(className,".R", sep=""),open="a")
		writeLines(src,con=outFile)
		close(inFile)
	}
	
	className<-paste("mtk",name,"DesignerResult", sep="")
	cat(file=outFile, "# Mexico Toolkit\n")
	cat(file=outFile, "# Author(s) :", authors,"\n")
	cat(file=outFile, "# Repository: https://mulcyber.toulouse.inra.fr/projects/baomexico/\n")
	cat(file=outFile, "# This file was generated automatically by the tool: mtk.AnalyserAddons(), which is built by Juhui WANG, INRA-MIA, Jouy.\n\n\n")
	if(!is.null(library))cat(file=outFile, "library(",library,")\n\n")
	
	cat(file=outFile, "#' A sub-class of the class \\code{\\linkS4class{mtkDesignerResult}} used to hold the results of the experiment design\n")
	if(!is.null(where))cat(file=outFile, "#' with the \"", name, "\" method defined in the \"", where, "\" file.\n", sep="")
	if(!is.null(library)) cat(file=outFile, "#' with the \"", name, "\"  method defined in the \"", library, "\" library.\n", sep="")
	if(!is.null(where))cat(file=outFile, "#' For more details, see the help of the function \"", main, "()\" defined in the \"", where, "\" file.\n", sep="")
	if(!is.null(library)) cat(file=outFile, "#' For more details, see the help of the function \"", main, "()\" defined in the  \"", library, "\" library .\n", sep="")
	
	cat(file=outFile, "#' @title The mtk", name, "DesignerResult class\n", sep="")
	cat(file=outFile, "#' @exportClass mtk",name,"DesignerResult\n\n", sep="")
	
	cat(file=outFile, "setClass(\"",className,"\",\n", sep="")
	cat(file=outFile, "	
					contains=c(\"mtkDesignerResult\")\n		)\n\n")
	
	cat(file=outFile, "#' The constructor.
#'  @param main a data frame to hold the main results produced by the Designer.
#'  @param information a named list to provide supplementary information about the sampling process and its results.
					\n")
	
	cat(file=outFile, "#' @return an object of class \\code{\\linkS4class{mtk",name,"DesignerResult}}\n",sep="")
	cat(file=outFile, "#' @examples mtk",className,"()\n",sep="")
	cat(file=outFile, "#' @export mtk",className,"\n",sep="")
	cat(file=outFile, "#' @title The constructor\n\n")
	
	cat(file=outFile, "mtk",name,"DesignerResult <- function(main, information=NULL) {\n",sep="")
	cat(file=outFile, "	res <- new(\"",className,"\", main=main, information=information)
					return(res)
					}\n", sep="")
	
	
	
	if(!is.null(summary)){
		cat(file=outFile, "#' Shows a summary of the results  produced by the Designer.\n")
		if(!is.null(where))cat(file=outFile, "#' For more details, see the help of the function \"", summary, "()\" defined in the \"", where, "\" file.\n", sep="")
		if(!is.null(library)) cat(file=outFile, "#' For more details, see the help of the function \"", summary, "()\" defined in the \"", library, "\" library.\n", sep="")
		
		cat(file=outFile, "#' @title The summary method\n")
		cat(file=outFile, "#' @param object an object of class \\code{\\linkS4class{mtk", name,"DesignerResult}}\n",sep="")
		cat(file=outFile, "
#' @return invisible()
#' @exportMethod summary\n\n")
		
	cat(file=outFile, "setMethod(f=\"summary\", \"", className, "\",\n", sep="")
	cat(file=outFile, "	definition=function(object,...){\n")
	
	cat(file=outFile, summary,"(list(main=object@main,information=object@information),...)\n",sep="")
	cat(file=outFile, "	})\n\n")
	}
	if(!is.null(plot)){
		cat(file=outFile, "#' Plots the results  produced by the Designer.\n")
		if(!is.null(where))cat(file=outFile, "#' For more details, see the help of the function \"", plot, "()\" defined in the \"", where, "\" file.\n", sep="")
		if(!is.null(library)) cat(file=outFile, "#' For more details, see the help of the function \"", plot, "()\" defined in the \"", library, "\" library.\n", sep="")
		
		cat(file=outFile, "#' @title The plot method\n")
		cat(file=outFile, "#' @param x an object of class \\code{\\linkS4class{mtk", name,"DesignerResult}}\n",sep="")
		cat(file=outFile, "
#' @return invisible()
#' @exportMethod plot\n\n")
	
		cat(file=outFile, "setMethod(f=\"plot\", \"", className, "\",\n", sep="")
		cat(file=outFile, "	definition=function(x,y,...){\n")
		cat(file=outFile, "if(!missing(y)) ", plot,"(list(main=x@main,information=x@information), y, ...)\n",sep="")
		cat(file=outFile, "else ", plot,"(list(main=x@main,information=x@information), ...)\n",sep="")
		cat(file=outFile, "	})\n")	
		}
	
	if(!is.null(print)){
		cat(file=outFile, "#' Prints the results  produced by the Designer.\n")
		if(!is.null(where))cat(file=outFile, "#' For more details, see the help of the function \"", print, "()\" defined in the \"", where, "\" file.\n", sep="")
		if(!is.null(library)) cat(file=outFile, "#' For more details, see the help of the function \"", print, "()\" defined in the \"", library, "\" library .\n", sep="")
		
		cat(file=outFile, "#' @title The print method\n")
		cat(file=outFile, "#' @param object an object of class \\code{\\linkS4class{mtk", name,"DesignerResult}}\n",sep="")
		cat(file=outFile, "
#' @return invisible()
#' @exportMethod print\n\n")
		
		cat(file=outFile, "setMethod(f=\"print\", \"", className, "\",\n", sep="")
		cat(file=outFile, "	definition=function(x, ...){\n")
		cat(file=outFile, print,"(list(main=x@main,information=x@information), ...)\n",sep="")
		cat(file=outFile, "	})\n")	
		}	
	
	close(outFile)
	}	


	

#'  A fonction used to extend the MTK package
#'  with existing models programmed as R functions. The native models must have the following syntax:
#'  main <- function(X, ...)
#'  where X is a data frame holding the experiment design used to run the simulation.
#' 
#'  Otherwise, the native model must return a named list with two elements: the element "main" holds the result of the simulation and
#'  the element "information" gives supplementary information about the simulation and its results.
#' 
#' 	Users can redefine the generic functions: summary (object, ...), plot <- function(x,y, ...), print(x, ...).
#' 
#' 
#'  For example, a model programmed by H. Monod is held in file "wwdm.R". The main function is named as wwdm.simule(). To convert the file "wwdm.R" into 
#'  mtk compliant S4 classes, we use the following call:
#' 	mtk.EvaluatorAddons(where="wwdm.R", authors="R. Monod, INRA-MIA Jouy", name="HMWWDM", main="wwdm.simule")
#'	This function generates a S4 class "mtkWWDMEvaluator.R" which can be directly integrated into the mtk package.

mtk.evaluatorAddons = function(where=NULL, library=NULL, authors=NULL, name=NULL, main=NULL,summary=NULL,plot=NULL, print=NULL)
{
	
	if(is.null(library)&&!file.exists(where))stop(paste("The file: ", where," or the library: ",library, "does not exist in this context!", sep=" "))
	className<-paste("mtk",name,"Evaluator", sep="")
	outFile<-file(paste(className,".R", sep=""),open="w")
	cat(file=outFile, "# Mexico Toolkit\n")
	cat(file=outFile, "# Author(s) :", authors,"\n")
	cat(file=outFile, "# Repository: https://mulcyber.toulouse.inra.fr/projects/baomexico/\n")
	cat(file=outFile, "# This file was generated automatically by the tool: mtk.EvaluatorAddons(), which is built by Juhui WANG, INRA-MIA, JOUY.\n\n\n")
	if(!is.null(library))cat(file=outFile, "library(",library,")\n\n")
	
	cat(file=outFile, "#' A sub-class of the class \\code{\\linkS4class{mtkEvaluator}} used to perform model simulation\n")
	if(!is.null(where))cat(file=outFile, "#' with the \"", name, "\" model defined in the \"", where, "\" file.\n", sep="")
	if(!is.null(library)) cat(file=outFile, "#' with the \"", name, "\"  model defined in the \"", library, "\" library.\n", sep="")
	if(!is.null(where))cat(file=outFile, "#' For more details, see the help of the function \"", main, "()\" defined in the \"", where, "\" file.\n", sep="")
	if(!is.null(library)) cat(file=outFile, "#' For more details, see the help of the function \"", main, "()\" defined in the  \"", library, "\" library .\n", sep="")
	
	cat(file=outFile, "#' @title The mtk", name, "Evaluator class\n", sep="")
	cat(file=outFile, "#' @exportClass mtk",name,"Evaluator\n\n", sep="")
	
	cat(file=outFile, "setClass(\"mtk",name,"Evaluator\",\n", sep="")
	cat(file=outFile, "		contains=c(\"mtkEvaluator\")\n")
	cat(file=outFile, "	)\n\n")
	
	cat(file=outFile, "#' The constructor.
#' @param mtkParameters a vector of [\\code{\\linkS4class{mtkParameter}}] representing the parameters necessary to run the Evaluator.
#' @param listParameters a named list defining the parameters necessary to run the Evaluator. It gives non object-oriented way to define the parameters.\n")
	
	cat(file=outFile, "#' @return an object of class \\code{\\linkS4class{mtk",name,"Evaluator}}\n",sep="")
	cat(file=outFile, "#' @examples mtk",name,"Evaluator()\n",sep="")
	cat(file=outFile, "#' @export mtk",name,"Evaluator\n",sep="")
	cat(file=outFile, "#' @title The constructor\n\n")
	
	cat(file=outFile, className,"<- function(mtkParameters=NULL, listParameters=NULL) {\n",sep="")
	cat(file=outFile, "	p<-mtkParameters
		if(!is.null(listParameters))
			p <- make.mtkParameterList(listParameters)
					\n")
	cat(file=outFile, "	res <- new(\"",className,"\",service=\"",name,"\", parameters=p)
		return(res)
		}\n", sep="")

	if(!is.null(where))cat(file=outFile, "\n\n#' Performs the simulation  with the model  \"",name,"\" defined in the \"",where,"\" file. \n", sep="")
	if(!is.null(library)) cat(file=outFile, "\n\n#' Performs the simulation  with the model  \"",name,"\" defined in the \"",library,"\" library.\n", sep="")
	cat(file=outFile, "#' @title The run method\n")
	cat(file=outFile, "#' @param this an object of class \\code{\\linkS4class{mtk", name,"Evaluator}}\n",sep="")
	cat(file=outFile, "#' @param context an object of class \\code{\\linkS4class{mtkExpWorkflow}}
#' @return invisible()
#' @exportMethod run
	\n")

	cat(file=outFile, "setMethod(f=\"run\", signature=c(this=\"", className, "\",\n", sep="")
	cat(file=outFile, "			context=\"mtkExpWorkflow\"),
		definition=function(this, context){
			if(this@state) return(invisible())
			nameThis<-deparse(substitute(this))
			
			X <- context@processesVector$design@result@main 
			parameters<-getParameters(this)
			##!!
			##!! Pre-processing the input data, the processing of the method to implement follows:
			##!!
	\n")
	
	cat(file=outFile, "	
					output<-eval(do.call(\"",main,"\",c(list(X=X), parameters)))\n",sep="")
		
	cat(file=outFile, "
			##!!
			##!!  post-processing the output of the method:
			##!!
	\n")
	cat(file=outFile, "		this@result <- mtk",name,"EvaluatorResult(main=output$main, information=output$information)\n",sep="")
	cat(file=outFile, "		this@state<-TRUE
			
			assign(nameThis, this, envir=parent.frame())
			return(invisible())
		})
	\n")
if(!is.null(where))
{
	cat(file=outFile, "				####################### \n")
	cat(file=outFile, "				## THE ORIGINAL CODE ## \n")
	cat(file=outFile, "				####################### \n\n\n")
	
	inFile<-file(where, open="r")
	src<-readLines(con=inFile, n=-1)
	close(outFile)
	outFile<-file(paste(className,".R", sep=""),open="a")
	writeLines(src,con=outFile)
	close(inFile)
}


className<-paste("mtk",name,"EvaluatorResult", sep="")
cat(file=outFile, "# Mexico Toolkit\n")
cat(file=outFile, "# Author(s) :", authors,"\n")
cat(file=outFile, "# Repository: https://mulcyber.toulouse.inra.fr/projects/baomexico/\n")
cat(file=outFile, "# This file was generated automatically by the tool: mtk.EvaluatorAddons() built by Juhui WANG, INRA-MIA, JOUY en Josas, FRANCE.\n\n\n")
if(!is.null(library))cat(file=outFile, "library(",library,")\n\n")

cat(file=outFile, "#' A sub-class of the class \\code{\\linkS4class{mtkEvaluatorResult}} used to hold the results of the model simulation \n")
if(!is.null(where))cat(file=outFile, "#' with the \"", name, "\" model defined in the \"", where, "\" file.\n", sep="")
if(!is.null(library)) cat(file=outFile, "#' with the \"", name, "\"  model defined in the \"", library, "\" library.\n", sep="")
if(!is.null(where))cat(file=outFile, "#' For more details, see the help of the function \"", main, "()\" defined in the \"", where, "\" file.\n", sep="")
if(!is.null(library)) cat(file=outFile, "#' For more details, see the help of the function \"", main, "()\" defined in the  \"", library, "\" library .\n", sep="")

cat(file=outFile, "#' @title The mtk", name, "EvaluatorResult class\n", sep="")
cat(file=outFile, "#' @exportClass mtk",name,"EvaluatorResult\n\n", sep="")

cat(file=outFile, "setClass(\"",className,"\",\n", sep="")
cat(file=outFile, "	
				contains=c(\"mtkEvaluatorResult\")\n		)\n\n")

cat(file=outFile, "#' The constructor.
#'  @param main a data-frame to hold the main results produced by the Evaluator.
#'  @param information a named list to provide supplementary information about the model simulation  and its results.
				\n")

cat(file=outFile, "#' @return an object of class \\code{\\linkS4class{mtk",name,"EvaluatorResult}}\n",sep="")
cat(file=outFile, "#' @examples mtk",className,"()\n",sep="")
cat(file=outFile, "#' @export mtk",className,"\n",sep="")
cat(file=outFile, "#' @title The constructor\n\n")

cat(file=outFile, "mtk",name,"EvaluatorResult <- function(main, information=NULL) {\n",sep="")
cat(file=outFile, "	res <- new(\"",className,"\", main=main, information=information)
				return(res)
				}\n", sep="")

	
	if(!is.null(summary)){
	cat(file=outFile, "#' Shows a summary of the results  produced by the Evaluator.\n")
	if(!is.null(where))cat(file=outFile, "#' For more details, see the help of the function \"", summary, "()\" defined in the \"", where, "\" file.\n", sep="")
	if(!is.null(library)) cat(file=outFile, "#' For more details, see the help of the function \"", summary, "()\" defined in the \"", library, "\" library.\n", sep="")
	
	cat(file=outFile, "#' @title The summary method\n")
	cat(file=outFile, "#' @param object an object of class \\code{\\linkS4class{mtk", name,"EvaluatorResult}}\n",sep="")
	cat(file=outFile, "#' @return invisible()
#' @exportMethod summary\n\n")

	cat(file=outFile, "setMethod(f=\"summary\", \"", className, "\",\n", sep="")
	cat(file=outFile, "	definition=function(object,...){\n")
			
	cat(file=outFile, summary,"(list(main=object@main,information=object@information),...)\n",sep="")
	cat(file=outFile, "	})\n\n")
	}
	if(!is.null(plot)){
	cat(file=outFile, "#' Plots the results  produced by the Evaluator.\n")
	if(!is.null(where))cat(file=outFile, "#' For more details, see the help of the function \"", plot, "()\" defined in the \"", where, "\" file.\n", sep="")
	if(!is.null(library)) cat(file=outFile, "#' For more details, see the help of the function \"", plot, "()\" defined in the \"", library, "\" library.\n", sep="")
	
	cat(file=outFile, "#' @title The plot method\n")
	cat(file=outFile, "#' @param x an object of class \\code{\\linkS4class{mtk", name,"EvaluatorResult}}\n",sep="")
	cat(file=outFile, "
#' @return invisible()
#' @exportMethod plot\n\n")

	cat(file=outFile, "setMethod(f=\"plot\", \"", className, "\",\n", sep="")
	cat(file=outFile, "	definition=function(x,y,...){\n")
	cat(file=outFile, "if(!missing(y)) ", plot,"(list(main=x@main,information=x@information),y, ...)\n",sep="")
	cat(file=outFile, "else ", plot,"(list(main=x@main,information=x@information), ...)\n",sep="")
	cat(file=outFile, "	})\n")	
	}
	if(!is.null(print)){
	cat(file=outFile, "#' Prints the results  produced by the Evaluator.\n")
	if(!is.null(where))cat(file=outFile, "#' For more details, see the help of the function \"", print, "()\" defined in the \"", where, "\" file.\n", sep="")
	if(!is.null(library)) cat(file=outFile, "#' For more details, see the help of the function \"", print, "()\" defined in the \"", library, "\" library.\n", sep="")
	
	cat(file=outFile, "#' @title The print method\n")
	cat(file=outFile, "#' @param object an object of class \\code{\\linkS4class{mtk", name,"EvaluatorResult}}\n",sep="")
	cat(file=outFile, "#' @return invisible()
#' @exportMethod print\n\n")
	
	cat(file=outFile, "setMethod(f=\"print\", \"", className, "\",\n", sep="")
	cat(file=outFile, "	definition=function(x, ...){\n")
	cat(file=outFile, print,"(list(main=x@main,information=x@information), ...)\n",sep="")
	cat(file=outFile, "	})\n")	
	}	
	
	close(outFile)
	
	}

#'  A fonction used to extend the MTK package
#'  with existing sensitivity analysis methods programmed as R functions. The Analyser must have the following syntax:
#'  main <- function(X,Y, ...)
#'  where the data frame "X" keeps the experiment design, and the data frame "Y" 
#'  holds the results of the model simulation.
#' 
#'  The function "main" returns the analyses results as a named list with two elements: the element "main" holds the result of the analysis and
#'  the element "information" gives supplementary information about the analysis process and its results.
#' 
#' 	Furthermore, users can redefine the generic functions: summary (object, ...), plot <- function(x,y, ...), print(x, ...).
#' 
#' 
#'  For example, the sensitivity analysis method "PLMM" programmed by R. Faivre is held in file "plmm-mtk.R". It contains three functions:
#'  the main function "plmm.mtk()", the summary function "summary.plmm()" and the plot function "plot.plmm()"
#'	To convert the file "plmm-mtk.R" into  mtk compliant classes, we use the following call:
#' 	mtk.AnalyserAddons(where="plmm-mtk.R", authors="R. Faivre, INRA-MIA Toulouse", name="PLMM", main="plmm.mtk", summary="summary.plmm", plot="plot.plmm")
#'	This function generates a mtk compliant S4 class named "mtkPLMMAnalyser.R" which can be directly integrated into the mtk package.

mtk.analyserAddons = function(where=NULL, library=NULL, authors=NULL, name=NULL, main=NULL,summary=NULL,plot=NULL, print=NULL)
{
	
	if(is.null(library)&&!file.exists(where))stop(paste("The file: ", where," or the library: ",library, "does not exist in this context!", sep=" "))
	className<-paste("mtk",name,"Analyser", sep="")
	outFile<-file(paste(className,".R", sep=""),open="w")
	cat(file=outFile, "# Mexico Toolkit\n")
	cat(file=outFile, "# Author(s) :", authors,"\n")
	cat(file=outFile, "# Repository: https://mulcyber.toulouse.inra.fr/projects/baomexico/\n")
	cat(file=outFile, "# This file was generated automatically by the tool: mtk.AnalyserAddons() built by Juhui WANG, INRA-MIA, JOUY, FRANCE.\n\n\n")
	if(!is.null(library))cat(file=outFile, "library(",library,")\n\n")
	
	cat(file=outFile, "#' A sub-class of the class \\code{\\linkS4class{mtkAnalyser}} used to perform the sensitivity  analysis\n")
	if(!is.null(where))cat(file=outFile, "#' with the \"", name, "\" method defined in the \"", where, "\" file.\n", sep="")
	if(!is.null(library)) cat(file=outFile, "#' with the \"", name, "\"  method defined in the \"", library, "\" library.\n", sep="")
	if(!is.null(where))cat(file=outFile, "#' For more details, see the help of the function \"", main, "()\" defined in the \"", where, "\" file.\n", sep="")
	if(!is.null(library)) cat(file=outFile, "#' For more details, see the help of the function \"", main, "()\" defined in the  \"", library, "\" library .\n", sep="")
	
	cat(file=outFile, "#' @title The mtk", name, "Analyser class\n", sep="")
	cat(file=outFile, "#' @exportClass mtk",name,"Analyser\n\n", sep="")
	
	cat(file=outFile, "setClass(\"mtk",name,"Analyser\",\n", sep="")
	cat(file=outFile, "		contains=c(\"mtkAnalyser\")\n")
	cat(file=outFile, "	)\n\n")
	
	cat(file=outFile, "#' The constructor.
#' @param mtkParameters a vector of [\\code{\\linkS4class{mtkParameter}}] representing the parameters necessary to run the Analyser.
#' @param listParameters a named list defining the parameters necessary to run the Analyser. It gives non object-oriented way to specify the parameters.\n")
	
	cat(file=outFile, "#' @return an object of class \\code{\\linkS4class{mtk",name,"Analyser}}\n",sep="")
	cat(file=outFile, "#' @examples mtk",name,"Analyser()\n",sep="")
	cat(file=outFile, "#' @export mtk",name,"Analyser\n",sep="")
	cat(file=outFile, "#' @title The constructor\n\n")
	
	cat(file=outFile, className,"<- function(mtkParameters=NULL, listParameters=NULL) {\n",sep="")
	cat(file=outFile, "	p<-mtkParameters
		if(!is.null(listParameters))
			p <- make.mtkParameterList(listParameters)
					\n")
	cat(file=outFile, "	res <- new(\"",className,"\",service=\"",name,"\", parameters=p)
		return(res)
		}\n", sep="")

	if(!is.null(where))cat(file=outFile, "\n\n#' Performs  sensitivity analysis  with the method  \"",name,"\" defined in the \"",where,"\" file. \n", sep="")
	if(!is.null(library)) cat(file=outFile, "\n\n#' Performs  sensitivity analysis  with the method  \"",name,"\" defined in the \"",library,"\" library.\n", sep="")
	cat(file=outFile, "#' @title The run method\n")
	cat(file=outFile, "#' @param this an object of class \\code{\\linkS4class{mtk", name,"Analyser}}\n",sep="")
	cat(file=outFile, "#' @param context an object of class \\code{\\linkS4class{mtkExpWorkflow}}
#' @return invisible()
#' @exportMethod run
	\n")

	cat(file=outFile, "setMethod(f=\"run\", signature=c(this=\"", className, "\",\n", sep="")
	cat(file=outFile, "			context=\"mtkExpWorkflow\"),
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
	\n")
	
	cat(file=outFile, "	
					analysisOutput<-eval(do.call(\"",main,"\",c(list(X=X,Y=Y), parameters)))\n",sep="")
		
	cat(file=outFile, "
			##!!
			##!!  post-processing the output of the method:
			##!!
	\n")
	cat(file=outFile, "		this@result <- mtk",name,"AnalyserResult(main=analysisOutput$main, information=analysisOutput$information)\n",sep="")
	cat(file=outFile, "		this@state<-TRUE
			
			assign(nameThis, this, envir=parent.frame())
			return(invisible())
		})
	\n")
if(!is.null(where))
{
	cat(file=outFile, "				####################### \n")
	cat(file=outFile, "				## THE ORIGINAL CODE ## \n")
	cat(file=outFile, "				####################### \n\n\n")
	
	inFile<-file(where, open="r")
	src<-readLines(con=inFile, n=-1)
	close(outFile)
	outFile<-file(paste(className,".R", sep=""),open="a")
	writeLines(src,con=outFile)
	close(inFile)
}


className<-paste("mtk",name,"AnalyserResult", sep="")
cat(file=outFile, "# Mexico Toolkit\n")
cat(file=outFile, "# Author(s) :", authors,"\n")
cat(file=outFile, "# Repository: https://mulcyber.toulouse.inra.fr/projects/baomexico/\n")
cat(file=outFile, "# This file was generated automatically by the tool: mtk.AnalyserAddons() built by Juhui WANG, INRA-MIA, JOUY, France.\n\n\n")
if(!is.null(library))cat(file=outFile, "library(",library,")\n\n")

cat(file=outFile, "#' A sub-class of the class \\code{\\linkS4class{mtkAnalyserResult}} used to hold the results of the sensitivity  analysis\n")
if(!is.null(where))cat(file=outFile, "#' with the \"", name, "\" method defined in the \"", where, "\" file.\n", sep="")
if(!is.null(library)) cat(file=outFile, "#' with the \"", name, "\"  method defined in the \"", library, "\" library.\n", sep="")
if(!is.null(where))cat(file=outFile, "#' For more details, see the help of the function \"", main, "()\" defined in the \"", where, "\" file.\n", sep="")
if(!is.null(library)) cat(file=outFile, "#' For more details, see the help of the function \"", main, "()\" defined in the  \"", library, "\" library .\n", sep="")

cat(file=outFile, "#' @title The mtk", name, "AnalyserResult class\n", sep="")
cat(file=outFile, "#' @exportClass mtk",name,"AnalyserResult\n\n", sep="")

cat(file=outFile, "setClass(\"",className,"\",\n", sep="")
cat(file=outFile, "	
				contains=c(\"mtkAnalyserResult\")\n		)\n\n")

cat(file=outFile, "#' The constructor.
#'  @param main a data-frame to hold the main results produced by the Analyser.
#'  @param information a named list to provide supplementary information about the analysis process and its results.
				\n")

cat(file=outFile, "#' @return an object of class \\code{\\linkS4class{mtk",name,"AnalyserResult}}\n",sep="")
cat(file=outFile, "#' @examples mtk",className,"()\n",sep="")
cat(file=outFile, "#' @export mtk",className,"\n",sep="")
cat(file=outFile, "#' @title The constructor\n\n")

cat(file=outFile, "mtk",name,"AnalyserResult <- function(main, information=NULL) {\n",sep="")
cat(file=outFile, "	res <- new(\"",className,"\", main=main, information=information)
				return(res)
				}\n", sep="")

	
	if(!is.null(summary)){
	cat(file=outFile, "#' Shows a summary of the results  produced by the Analyser.\n")
	if(!is.null(where))cat(file=outFile, "#' For more details, see the help of the function \"", summary, "()\" defined in the \"", where, "\" file.\n", sep="")
	if(!is.null(library)) cat(file=outFile, "#' For more details, see the help of the function \"", summary, "()\" defined in the \"", library, "\" library.\n", sep="")
	
	cat(file=outFile, "#' @title The summary method\n")
	cat(file=outFile, "#' @param object an object of class \\code{\\linkS4class{mtk", name,"AnalyserResult}}\n",sep="")
	cat(file=outFile, "#' @return invisible()
#' @exportMethod summary\n\n")

	cat(file=outFile, "setMethod(f=\"summary\", \"", className, "\",\n", sep="")
	cat(file=outFile, "	definition=function(object,...){\n")
			
	cat(file=outFile, summary,"(list(main=object@main,information=object@information),...)\n",sep="")
	cat(file=outFile, "	})\n\n")
	}
	if(!is.null(plot)){
	cat(file=outFile, "#' Plots the results  produced by the Analyser.\n")
	if(!is.null(where))cat(file=outFile, "#' For more details, see the help of the function \"", plot, "()\" defined in the \"", where, "\" file.\n", sep="")
	if(!is.null(library)) cat(file=outFile, "#' For more details, see the help of the function \"", plot, "()\" defined in the \"", library, "\" library.\n", sep="")
	
	cat(file=outFile, "#' @title The plot method\n")
	cat(file=outFile, "#' @param x an object of class \\code{\\linkS4class{mtk", name,"AnalyserResult}}\n",sep="")
	cat(file=outFile, "
#' @return invisible()
#' @exportMethod plot\n\n")

	cat(file=outFile, "setMethod(f=\"plot\", \"", className, "\",\n", sep="")
	cat(file=outFile, "	definition=function(x,y,...){\n")
	cat(file=outFile, "if(!missing(y)) ", plot,"(list(main=x@main,information=x@information),y, ...)\n",sep="")
	cat(file=outFile, "else ", plot,"(list(main=x@main,information=x@information), ...)\n",sep="")
	cat(file=outFile, "	})\n")	
	}
	if(!is.null(print)){
	cat(file=outFile, "#' Prints the results  produced by the Analyser.\n")
	if(!is.null(where))cat(file=outFile, "#' For more details, see the help of the function \"", print, "()\" defined in the \"", where, "\" file.\n", sep="")
	if(!is.null(library)) cat(file=outFile, "#' For more details, see the help of the function \"", print, "()\" defined in the \"", library, "\" library.\n", sep="")
	
	cat(file=outFile, "#' @title The print method\n")
	cat(file=outFile, "#' @param object an object of class \\code{\\linkS4class{mtk", name,"AnalyserResult}}\n",sep="")
	cat(file=outFile, "#' @return invisible()
#' @exportMethod print\n\n")
	
	cat(file=outFile, "setMethod(f=\"print\", \"", className, "\",\n", sep="")
	cat(file=outFile, "	definition=function(x, ...){\n")
	cat(file=outFile, print,"(list(main=x@main,information=x@information), ...)\n",sep="")
	cat(file=outFile, "	})\n")	
	}	
	
	close(outFile)
	
	}