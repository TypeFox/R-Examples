# Mexico Toolkit
# Author(s) : H. Monod, INRA-MIA-Jouy, 78352, Jouy en Josas, France 
# Repository: https://mulcyber.toulouse.inra.fr/projects/baomexico/
# This file was generated automatically by the tool: mtk.designerAddons() built by Juhui WANG, INRA-JOUY, FRANCE.


#' A sub-class of the class \code{\linkS4class{mtkDesigner}} used to generate the experiment design
#' with the "Sobol" method defined in the "sobolDesigner.R" file.
#' For more details, see the help of the function "Designer.sobol()" defined in the "sobolDesigner.R" file.
#' @title The mtkSobolDesigner class
#' @exportClass mtkSobolDesigner

setClass("mtkSobolDesigner",
		contains=c("mtkDesigner")
	)

#' The constructor.
#' @param mtkParameters a vector of [\code{\linkS4class{mtkParameter}}] representing the parameters necessary to run the Designer.
#' @param listParameters a named list defining the parameters necessary to run the Designer. It gives non object-oriented way to define the value of the parameters.
#' @return an object of class \code{\linkS4class{mtkSobolDesigner}}
#' @examples mtkSobolDesigner()
#' @export mtkSobolDesigner
#' @title The constructor

mtkSobolDesigner <- function(mtkParameters=NULL, listParameters=NULL) {
	p<-mtkParameters
	if(!is.null(listParameters))
					p <- make.mtkParameterList(listParameters)
					
	res <- new("mtkSobolDesigner",service="Sobol", parameters=p)
	return(res)
					}


#' Generates the experiment design with the method  "Sobol" defined in the "sobolDesigner.R" file
#' @title The run method
#' @param this an object of class \code{\linkS4class{mtkSobolDesigner}}
#' @param context an object of class \code{\linkS4class{mtkExpWorkflow}}
#' @return invisible()
#' @exportMethod run
					
setMethod(f="run", signature=c(this="mtkSobolDesigner",
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
					sortie<-eval(do.call("Designer.sobol",c(arg, parameters)))

					##!!
					##!!  post-processing the output of the method:
					##!!
	
		this@result <- mtkSobolDesignerResult(main=sortie$main, information=sortie$information)
		this@state<-TRUE
					
					assign(nameThis, this, envir=parent.frame())
					return(invisible())
					})
					
				####################### 
				## THE ORIGINAL CODE ## 
				####################### 


#' A Designer using the Sobol'-Saltelli method of sensitivity analysis
#' @param factors a vector of factors' names or an integer giving the number of factors
#' @param distribNames a string or a character vector giving the names of the factors uncertainty distributions
#' @param distribParameters a list of lists of distribution parameters
#' @param N size of the basic samples; the final sample size will be N*(k+2)
#' @param sampling name of the basic sampling method to be used (presently, "MC" for Monte Carlo, "LHS" for Latin hypercube)
#' @param shrink a scalar or a vector of scalars between 0 and 1, specifying shrinkage to be used on the probabilities before calculating the quantiles. See the details.
#' @return a data.frame of N(k+2) rows and k columns
#' @note based on the sobol function of the sensitivity package
Designer.sobol <- function(factors, distribNames, distribParameters, N, sampling="MC", shrink=1,...){
  ## number of factors
  nbf <- length(factors)
  ## in case factors is an integer
  if((nbf==1) & is.numeric(factors)){
    nbf <- factors
    factors <- paste("X", seq(nbf), sep="")
  }
  ## distributions
  if(missing(distribNames)){
    distribNames <- rep("unif", nbf)
    distribParameters <- rep(list(list(min=0,max=1)), nbf)
  }
  if(length(distribNames) == 1){
    distribNames <- rep(distribNames, nbf)
    distribParameters <- rep(list(distribParameters), nbf)
  }
  
  ## Monte Carlo sampling
  if (sampling == "MC"){
    A <- data.frame(matrix(runif(N*nbf), N, nbf))
    B <- data.frame(matrix(runif(N*nbf), N, nbf))
  }
  ## Latin hypercube sampling
  else if (sampling == "LHS"){
    A <- data.frame(lhs003(nbf, N))
    B <- data.frame(lhs003(nbf, N))
  }
  ## 
  else (stop("Sampling method not included"))

  ## construction of the whole sample, based on the sensitivity package code
  ##X <- rbind(A, B)
  ##for (i in 1:k) {
  ##  Xb <- A
  ##  Xb[, i] <- B[, i]
  ##  X <- rbind(X, Xb)
  ##}

  ## Calculations: call to the morris function of the
  ## sensitivity library
  information <- list(X1=A, X2=B)
  sobol.sortie <- eval( do.call("sobol2002", information) )

  ## Post-processing, quantiles
  for(i in seq(nbf)){
    sobol.sortie$X[,i] <- Quantiles(sobol.sortie$X[,i],
				distribNames[i],
				distribParameters[[i]],
				shrink[i])
  }
  
  ## Output
  information <- sobol.sortie[c("model", "X1", "X2", "call")]

  information$SamplingMethod <- "sobol2002"
  information$InitialSampling <- sampling
  
  design <- as.data.frame(sobol.sortie$X)
  colnames(design) <- factors
  
  resultat <- list(main=design, information=information)
  
  return(resultat)
}

#' A simple function to generate Latin hypercubes
#' @param k number of factors
#' @param N size of the LHS
#' @return a data.frame
lhs003 <- function(k, N){
  X <- matrix(((0:(N-1)) + runif(N*k))/N, ncol=k, nrow=N)
  for(j in 1:k){
    X[,j] <- sample(X[,j],size=N)
  }
  as.data.frame(X)
}
# Mexico Toolkit
# Author(s) : H. Monod, INRA-MIA-Jouy, 78352, Jouy en Josas, France 
# Repository: https://mulcyber.toulouse.inra.fr/projects/baomexico/
# This file was generated automatically by the tool: mtk.AnalyserAddons(), which is built by Juhui WANG, INRA-MIA, Jouy.


#' A sub-class of the class \code{\linkS4class{mtkDesignerResult}} used to hold the results of the experiment design
#' with the "Sobol" method defined in the "sobolDesigner.R" file.
#' For more details, see the help of the function "Designer.sobol()" defined in the "sobolDesigner.R" file.
#' @title The mtkSobolDesignerResult class
#' @exportClass mtkSobolDesignerResult

setClass("mtkSobolDesignerResult",
	
					contains=c("mtkDesignerResult")
		)

#' The constructor.
#'  @param main a data frame to hold the main results produced by the Designer.
#'  @param information a named list to provide supplementary information about the sampling process and its results.
					
#' @return an object of class \code{\linkS4class{mtkSobolDesignerResult}}
#' @examples mtkmtkSobolDesignerResult()
#' @export mtkmtkSobolDesignerResult
#' @title The constructor

mtkSobolDesignerResult <- function(main=NULL, information=NULL) {
	res <- new("mtkSobolDesignerResult", main=main, information=information)
					return(res)
					}
