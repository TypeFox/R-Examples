##' MvBinary a package for Multivariate Binary data
##'
##' MvBinary is a tool for fitting the distribution of correlated multivariate binary data.
##'
##' \tabular{ll}{
##'   Package: \tab MvBinary\cr 
##'   Type: \tab Package\cr 
##'   Version: \tab 1.0.0\cr
##'   Date: \tab 2015-11-03\cr 
##'   License: \tab GPL-2\cr 
##'   LazyLoad: \tab yes\cr
##' }
##'
##'
##' @name MvBinary-package
##' @aliases MvBinary
##' @rdname MvBinary-package
##' @docType package
##' @keywords package
##' @import methods
##' @import mgcv
##' @import parallel
##' @import graphics
##' @export MvBinaryEstim
##' @export MvBinaryProbaPost
##' @export ComputeEmpiricCramer
##' @export ComputeMvBinaryCramer
##' @exportMethod print
##' @exportMethod summary
##' @exportClass MvBinaryResult
##' @importFrom stats as.dist cutree hclust optimize runif
##'
##' @author
##' Author: Marbac M., and Sedki S.
##'
##' @references Matthieu Marbac, Mohammed Sedki (2015). A Family of Blockwise One-Factor Distributions for Modelling High-Dimensional Binary Data. arXiv:1511.01343
##'
##' @examples
##' # Package loading
##' rm(list=ls())
##' require(MvBinary)
##' 
##' # Data loading
##' data(MvBinaryExample)
##' 
##' # Parameter estimation by the HAC-based algorithm on 2 cores
##' # where the EM algorithms are initialized 10 times
##' res.CAH <- MvBinaryEstim(MvBinaryExample, 2, nbinit.EM = 10)
##' 
##' # Summary of the estimated model
##' summary(res.CAH) 
##' 
##' # Print the parameters of the estimated model
##' print(res.CAH) 
NULL

##' Real binary data: Plants
##' 
##' The file plants.rda describes 35583 plants by indicating if they occur (1) or not (2) in 69 states of the Norht America.
##'
##' This data set been extracted from the USA plants database, July 29, 2015.
##'
##' @format A matrix with 35583 observations on the 69 variables.
##'
##' 
##'
##'
##' @name plants
##' @docType data
##' @keywords datasets
##' 
##' @examples
##'   data(plants)
NULL


##' Simulated binary data: MvBinaryExample
##' 
##' The file MvBinaryExample.rda describes 400 individuals by 6 binary variables.
##'
##' This data set has been simulated from the MvBinary model. The first three variables are dependent. The last three variables are dependent.
##'
##' @format A matrix with 400 observations on the 6 variables.
##'
##' 
##'
##'
##' @name MvBinaryExample
##' @docType data
##' @keywords datasets
##' 
##' @examples
##'   data(MvBinaryExample)
NULL




###################################################################################
##' Create an instance of the [\code{\linkS4class{MvBinaryResult}}] class
##'
##' This function performs the model selection and the parameter inference.
##' 
##' @param x matrix of the binary observation.
##' @param nbcores number of cores used for the model selection (only for Linux). Default is 1.
##' @param algorithm algorithm used for the model selection ("HAC": deterministic algorithm based on the HAC of the variables, "MH": stochastic algorithm for optimizing the BIC criterion, "List": comparison of the models provided by the users). Default is "HAC".
##' @param modelslist list of models provided by the user (only used when algorithm="List"). Default is NULL
##' @param tol.EM stopping criterion for the EM algorithm. Default is 0.01
##' @param nbinit.EM number of random initializations for the EM algorithm. Default is 40.
##' @param nbiter.MH number of successive iterations without finding a model having a better BIC criterion which involves the stopping of the Metropolis-Hastings algorithm (only used when algorithm="MH"). Default is 50.
##' @param nbchains.MH number of radom initializations for the stochastic algorithm (only used when algorithm="MH"). Default is 10.
##'
##' @examples
##' # Data loading
##' data(MvBinaryExample)
##' 
##' # Parameter estimation by the HAC-based algorithm on 2 cores
##' # where the EM algorithms are initialized 10 times
##' res.CAH <- MvBinaryEstim(MvBinaryExample, 2, nbinit.EM = 10)
##' 
##' # Parameter estimation for two competing models
##' res.CAH <- MvBinaryEstim(MvBinaryExample, algorithm="List",
##'  modelslist=list(c(1,1,2,2,3,4), c(1,1,1,2,2,2)), nbinit.EM = 10)
##' 
##' # Summary of the estimated model
##' summary(res.CAH) 
##' 
##' # Print the parameters of the estimated model
##' print(res.CAH) 
##'
##' @return Returns an instance of the [\code{\linkS4class{MvBinaryResult}}] class. 
##' @export
##'
##'
MvBinaryEstim <- function(x, nbcores=1, algorithm="HAC", modelslist=NULL, tol.EM=0.01, nbinit.EM=40, nbiter.MH=50, nbchains.MH=10){
  if ( (is.matrix(x)==FALSE) || any((x==0) + (x==1) == 0) )
    stop("The input parameter x must be a binary matrix")
  
  if (is.null(colnames(x))) 
    colnames(x) <- paste("x",1:ncol(x), sep="")  
  
  if ((is.numeric(nbcores)==FALSE) || (length(nbcores)!=1) || (nbcores!=ceiling(nbcores)) )
    stop("The input parameter nbcores must be an integer")
  
  if (algorithm %in% c("HAC", "MH", "List") == FALSE)
    stop("The input parameter algorithm must take one of these values HAC, MH, List")
  
  if (is.null(modelslist)==FALSE){
    if (is.list(modelslist)){
      for (j in 1:length(modelslist)){
        if (length(modelslist[[j]])!=ncol(x))       stop("The input parameter modelslist must be a list where each element gives the partition of the variables by a vector of size d")
      }
    }else{
      stop("The input parameter modelslist must be a list where each element gives the partition of the variables by a vector of size d")
    }
  }
  
  if ((is.numeric(tol.EM)==FALSE) || (length(tol.EM)!=1))
    stop("The input parameter tol.EM must be a numeric value")
  
  if ((is.numeric(nbinit.EM)==FALSE) || (length(nbinit.EM)!=1) || (nbinit.EM!=ceiling(nbinit.EM)) )
    stop("The input parameter nbinit.EM must be an integer")
  
  if (algorithm=="MH")
    output <- MvBinaryEstimMH(x, nbcores, tol.EM, nbinit.EM, nbiter.MH, nbchains.MH)
  else
    output <- MvBinaryEstimCAH(x, nbcores, tol.EM, nbinit.EM, modelslist)
  return(output)
}

MvBinaryEstimCAH <- function(x, nbcores=1, tol=0.01, nbinit.EM=40, modelslist=NULL){
  alpha <- colMeans(x)
  if (is.null(modelslist)){
    # Computation of the Cramer's V 
    VcramerEmpiric <- ComputeEmpiricCramer(x)
    tree <- hclust(as.dist(1-VcramerEmpiric), method="ward")
    models <- list(); for (k in 1:ncol(x)) models[[k]] <- cutree(tree, k)
  }else{
      models <- modelslist
  }
  # Inference for the competiting models
  nb.cpus <- min(detectCores(all.tests = FALSE, logical = FALSE), nbcores)
  if ((nbcores>1)&&(Sys.info()["sysname"] != "Windows")){
    reference <- mclapply(X = models,
                          FUN = XEMmodel,
                          dataset=x,
                          alpha=alpha,
                          tol=tol,
                          nbinit.EM=nbinit.EM,
                          mc.cores = nb.cpus, mc.preschedule = TRUE, mc.cleanup = TRUE)
  }else{
    reference <- list(); for (loc in 1:length(models)) reference[[loc]] <- XEMmodel(x, alpha, tol, nbinit.EM, models[[loc]])
  }
  # Design outputs
  allBIC <- rep(NA, length(reference))
  allModels <- matrix(NA, length(reference), ncol(x))
  for (loc in 1:length(reference)){
    allBIC[loc] <- reference[[loc]]$bic
    allModels[loc,] <- models[[loc]]
  }
  Best <- reference[[which.max(allBIC)]]
  names(reference[[which.max(allBIC)]]$epsilon) <-  names(reference[[which.max(allBIC)]]$delta) <- colnames(x)
  return( new("MvBinaryResult", 
              alpha=alpha, 
              epsilon=reference[[which.max(allBIC)]]$epsilon, 
              delta=reference[[which.max(allBIC)]]$delta, 
              blocks=reference[[which.max(allBIC)]]$model,
              nbparam=reference[[which.max(allBIC)]]$nbparam,
              loglike=reference[[which.max(allBIC)]]$loglike,
              bic=reference[[which.max(allBIC)]]$bic)
  )
}

MvBinaryEstimMH <- function(x, nbcores=1, tol.EM=0.01, nbinit.EM=40, nbiter.MH=50, nbchains.MH=10){
  alpha <- colMeans(x)
  # Inference for the competiting models
  nb.cpus <- min(detectCores(all.tests = FALSE, logical = FALSE), nbcores, nbchains.MH)
  if ((nbcores>1)&&(Sys.info()["sysname"] != "Windows")){
    reference <- mclapply(X = as.list(rep(nbiter.MH, nbchains.MH)),
                          FUN = OneMH,
                          x=x,
                          alpha=alpha,
                          tol=tol.EM,
                          nbinit=nbinit.EM,
                          mc.cores = nb.cpus,
                          mc.preschedule = TRUE,
                          mc.cleanup = TRUE)
  }else{
    reference <- list(); for (loc in 1:nbchains.MH) reference[[loc]] <- OneMH(x, alpha, tol.EM, nbinit.EM, nbiter.MH)
  }
  # Design outputs
  allBIC <- rep(NA, length(reference))
  for (loc in 1:length(reference))    allBIC[loc] <- reference[[loc]]$bic
  Best <- XEMmodel(x, alpha, tol.EM, nbinit.EM, reference[[which.max(allBIC)]]$blocks)
  names(Best$epsilon) <-  names(Best$delta) <- colnames(x)
  return( new("MvBinaryResult", 
              alpha=alpha, 
              epsilon=Best$epsilon, 
              delta=Best$delta, 
              blocks=Best$model,
              nbparam=Best$nbparam,
              loglike=Best$loglike,
              bic=Best$bic)
  )
}