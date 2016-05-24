# Mexico Toolkit
# Author(s) : H. Monod, INRA-MIA-Jouy, 78352, Jouy en Josas, France 
# Repository: https://mulcyber.toulouse.inra.fr/projects/baomexico/
# This file was generated automatically by the tool: mtk.designerAddons() built by Juhui WANG, INRA-JOUY, FRANCE.


#' A sub-class of the class \code{\linkS4class{mtkDesigner}} used to generate the experiment design
#' with the "Morris" method defined in the "morrisSampler.R" file.
#' For more details, see the help of the function "sampler.morris()" defined in the "morrisSampler.R" file.
#' @title The mtkMorrisDesigner class
#' @exportClass mtkMorrisDesigner

setClass("mtkMorrisDesigner",
		contains=c("mtkDesigner")
	)

#' The constructor.
#' @param mtkParameters a vector of [\code{\linkS4class{mtkParameter}}] representing the parameters necessary to run the Designer.
#' @param listParameters a named list defining the parameters necessary to run the Designer. It gives non object-oriented way to define the value of the parameters.
#' @return an object of class \code{\linkS4class{mtkMorrisDesigner}}
#' @examples mtkMorrisDesigner()
#' @export mtkMorrisDesigner
#' @title The constructor

mtkMorrisDesigner <- function(mtkParameters=NULL, listParameters=NULL) {
	p<-mtkParameters
	if(!is.null(listParameters))
					p <- make.mtkParameterList(listParameters)
					
	res <- new("mtkMorrisDesigner",service="Morris", parameters=p)
	return(res)
					}


#' Generates the experiment design with the method  "Morris" defined in the "morrisSampler.R" file
#' @title The run method
#' @param this an object of class \code{\linkS4class{mtkMorrisDesigner}}
#' @param context an object of class \code{\linkS4class{mtkExpWorkflow}}
#' @return invisible()
#' @exportMethod run
					
setMethod(f="run", signature=c(this="mtkMorrisDesigner",
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
					sortie<-eval(do.call("sampler.morris",c(arg, parameters)))

					##!!
					##!!  post-processing the output of the method:
					##!!
	
		this@result <- mtkMorrisDesignerResult(main=sortie$main, information=sortie$information)
		this@state<-TRUE
					
					assign(nameThis, this, envir=parent.frame())
					return(invisible())
					})
					
				####################### 
				## THE ORIGINAL CODE ## 
				####################### 


#' A sampler built with the Morris method and its implementation in the sensitivity library, generalised to any set of uncertainty distributions
#' @param factors vector of factors names
#' @param distribNames a string or a character vector giving the names of the factors uncertainty distributions
#' @param distribParameters a list of lists of distribution parameters
#' @param type see the help on function morris in the sensitivity library
#' @param levels see the help on function morris in the sensitivity library
#' @param grid.jump see the help on function morris in the sensitivity library
#' @param r see the help on function morris in the sensitivity library
#' @param scale.factor see the help on function morris in the sensitivity library
#' @param scale see the help on function morris in the sensitivity library
#' @param identify see the help on function morris in the sensitivity library
#' @param shrink a scalar or a vector of scalars between 0 and 1, specifying shrinkage to be used on the probabilities before calculating the quantiles. See the details.
#' @value a list consisting of a design data.frame and a list of information parameters
#' @note The sampling design is calculated in two steps: first, an initial design assuming uniform 0-1 uncertainty distributions is calculated by the morris function of the sensitivity library. Then the final design is adapted to the factors uncertainty distributions by calculating quantiles based on the probabilities given by the initial design. When the uncertainty distributions are not bounded, the quantiles associated with 0 and 1 are infinite. To avoid that, it is necessary to set the shrink argument to a value strictly smaller than 1. 

sampler.morris <- function(factors, distribNames, distribParameters, type, levels, grid.jump, r, scale.factor=TRUE, scale=TRUE, identify, shrink=1, ...){
  ## number of factors
  nbf <- length(factors)
  ## in case factors is an integer
  if((nbf==1) & is.numeric(factors)){
    nbf <- factors
    factors <- paste("X", seq(nbf), sep="")
  }
  ##
  if(missing(distribNames)){
    distribNames <- rep("unif", nbf)
    distribParameters <- list(min=0,max=1)
  }
  if(length(distribNames) == 1){
    distribNames <- rep(distribNames, nbf)
    distribParameters <- rep(list(distribParameters), nbf)
 }
  
  ## Pre-processing on the process parameters to get the arguments as required by the morris function
  ## Case 1: a classical morris design is required

   if(type=="oat"){
    if(is.null(grid.jump)){
      grid.jump <- floor(levels/2)  ## default suggested by Saltelli et al
    }
    morrisDesignParameters <- list(type="oat",
                                   levels=levels,
                                   grid.jump=grid.jump
                                   )
  }
  ## Case 2: a "Pujol" simplex design is required
  if(type=="simplex"){
    morrisDesignParameters <- list(type="simplex",
                                   scale.factor=scale.factor
                                   )
  }
  ## End of the pre-processing
  morrisProcessParameters <- list(r=r,
                                  design=morrisDesignParameters, scale=scale, ...)
  
  ## Calculations: call to the morris function of the
  ## sensitivity library
  information <- c(list(factors=factors), morrisProcessParameters)
  morris.sortie <- eval( do.call("morris", information) )
  
  ## quantile calculations
  binf <- rep(NA,length=nbf)
  bsup <- rep(NA,length=nbf)
	
  for(i in seq(nbf)){
    morris.sortie$X[,i] <- Quantiles(morris.sortie$X[,i],
				distribNames[i],
				distribParameters[[i]],
				shrink[i])
    binf[i] <- Quantiles(0, distribNames[i], distribParameters[[i]], shrink[i])
    bsup[i] <- Quantiles(1, distribNames[i], distribParameters[[i]], shrink[i])
  }
  
  ## Output
  information <- morris.sortie[c("model", "factors", "r",
                                 "design", "binf", "bsup", 
                                 "scale", "call")]
  names(information)[3] <- "r.design" # to avoid conflicts in the Analyser
  information$binf <- binf
  information$bsup <- bsup

  information$SamplingMethod <- "morris"
  
  design <- as.data.frame(morris.sortie$X)
  
  resultat <- list(main=design, information=information)
  
  return(resultat)
}
# Mexico Toolkit
# Author(s) : H. Monod, INRA-MIA-Jouy, 78352, Jouy en Josas, France 
# Repository: https://mulcyber.toulouse.inra.fr/projects/baomexico/
# This file was generated automatically by the tool: mtk.AnalyserAddons(), which is built by Juhui WANG, INRA-MIA, Jouy.


#' A sub-class of the class \code{\linkS4class{mtkDesignerResult}} used to hold the results of the experiment design
#' with the "Morris" method defined in the "morrisSampler.R" file.
#' For more details, see the help of the function "sampler.morris()" defined in the "morrisSampler.R" file.
#' @title The mtkMorrisDesignerResult class
#' @exportClass mtkMorrisDesignerResult

setClass("mtkMorrisDesignerResult",
	
					contains=c("mtkDesignerResult")
		)

#' The constructor.
#'  @param main a data frame to hold the main results produced by the Designer.
#'  @param information a named list to provide supplementary information about the sampling process and its results.
					
#' @return an object of class \code{\linkS4class{mtkMorrisDesignerResult}}
#' @examples mtkmtkMorrisDesignerResult()
#' @export mtkmtkMorrisDesignerResult
#' @title The constructor

mtkMorrisDesignerResult <- function(main, information=NULL) {
	res <- new("mtkMorrisDesignerResult", main=main, information=information)
					return(res)
					}
