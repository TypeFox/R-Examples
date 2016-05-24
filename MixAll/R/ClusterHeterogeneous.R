#-----------------------------------------------------------------------
#     Copyright (C) 2012-2014  Serge Iovleff, University Lille 1, Inria
#
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as
#    published by the Free Software Foundation; either version 2 of the
#    License, or (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public
#    License along with this program; if not, write to the
#    Free Software Foundation, Inc.,
#    59 Temple Place,
#    Suite 330,
#    Boston, MA 02111-1307
#    USA
#
#    Contact : S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
#
#-----------------------------------------------------------------------
#' @include ClusterModelNames.R IClusterModel.R
NULL

#-----------------------------------------------------------------------
#' Create an instance of the [\code{\linkS4class{ClusterHeterogeneous}}] class
#'
#' This function computes the optimal heterogeneous mixture model according
#' to the \code{criterion} among the number of clusters given in
#' \code{nbCluster} using the strategy specified in [\code{strategy}].
#'
#' @param data [\code{list}] containing the data sets (matrices and/or data.frames).
#' If the data sets contains NA values, they will be estimated during the
#' estimation process.
#' @param modelNames [\code{vector}] of character and of same length than data
#' containing the model names to fit to each data set.
#' @param nbCluster  [\code{\link{vector}}] listing the number of clusters to test.
#' @param strategy a [\code{\linkS4class{ClusterStrategy}}] object containing
#' the strategy to run. clusterStrategy() method by default.
#' @param criterion character defining the criterion to select the best model.
#' The best model is the one with the lowest criterion value.
#' Possible values: "BIC", "AIC", "ICL". Default is "ICL".
#' @param nbCore integer defining the number of processors to use (default is 1, 0 for all).
#'
#' @examples
#' ## A quantitative example with the heart disease data set
#'  data(HeartDisease.cat)
#'  data(HeartDisease.cont)
#' ## with default values
#' ldata = list(HeartDisease.cat,HeartDisease.cont);
#' lnames = c("categorical_pk_pjk","gaussian_pk_sjk")
#' model <- clusterHeterogeneous(ldata, lnames, nbCluster=2:5, strategy = clusterFastStrategy())
#'
#' ## get summary
#' summary(model)
#'
#' ## print model
#' \dontrun{
#' print(model)
#' }
#' ## get estimated missing values
#' missingValues(model)
#'
#' ## use graphics functions
#' \dontrun{
#' plot(model)
#' }
#'
#' @return An instance of the [\code{\linkS4class{ClusterHeterogeneous}}] class.
#' @author Serge Iovleff
#' @export
#'
clusterHeterogeneous <- function( data, modelNames, nbCluster=2
                                , strategy=clusterFastStrategy()
                                , criterion="ICL"
                                , nbCore = 1)
{
  # check nbCluster
  nbClusterModel = length(nbCluster);
  nbClusterMin   = min(nbCluster);
  nbClusterMax   = max(nbCluster);
  if (nbClusterMin < 1) { stop("The number of clusters must be greater or equal to 1")}
  # check criterion
  if(sum(criterion %in% c("BIC","AIC", "ICL")) != 1)
  { stop("criterion is not valid. See ?clusterHeterogeneous for the list of valid criterion")}
  # check strategy
  if(class(strategy)[1] != "ClusterStrategy")
  {stop("strategy is not a Cluster Stategy class (must be an instance of the class ClusterStrategy).")}
  validObject(strategy);

  # check data and modelNames
  if (!is.list(data)) { stop("data must be a list");}
  if (!is.vector(modelNames)) { stop("modelNames must be a vector of character");}
  if (length(data) != length(modelNames)) { stop("data and modelNames must be equal lengths");}
  # create list of component
  ldata <- vector("list", length(data));
  for (i in 1:length(data))
  {
    if(clusterValidCategoricalNames(modelNames[i]) )
    { ldata[[i]] <- new("ClusterCategoricalComponent", data[[i]], nbClusterMin, modelNames[i]);}
    else
    {
      if( clusterValidGammaNames(modelNames[i]) )
      { ldata[[i]] <- new("ClusterGammaComponent", data[[i]], nbClusterMin, modelNames[i]);}
      else
      {
        if( clusterValidDiagGaussianNames(modelNames[i]) )
        { ldata[[i]] <- new("ClusterDiagGaussianComponent", data[[i]], nbClusterMin, modelNames[i]);}
        else
        stop("invalid model name");
      }
    }
  }
  # Create model
  model = new("ClusterHeterogeneous", ldata)
  model@strategy = strategy;
  # start estimation of the models
  resFlag  <- FALSE;
  if (length(nbCluster) >0)
  {
   resFlag = .Call("clusterMixtureHeterogene", model, nbCluster, strategy, criterion, nbCore, PACKAGE="MixAll");
  }
  # set names
  if (resFlag != TRUE) {cat("WARNING: An error occurs during the clustering process");}
  for (i in 1:length(data))
  {
    if(clusterValidCategoricalNames(modelNames[i]))
    { dim(model@ldata[[i]]@plkj) <- c(model@ldata[[i]]@nbModalities, model@nbCluster, ncol(model@ldata[[i]]@data))}
  }
  model
}

#-----------------------------------------------------------------------
#' Definition of the [\code{\linkS4class{ClusterHeterogeneous}}] class
#'
#' This class defines an Heterogeneous mixture Model.
#'
#' This class inherits from the [\code{\linkS4class{IClusterModelBase}}] class.
#' An heterogeneous model is a mixture model of the form:
#' \deqn{
#' f({{x}}_i=({{x}}_{1i}, {{x}}_{2i},\ldots {{x}}_{Li})|\theta)
#' = \sum_{k=1}^K p_k \prod_{l=1}^L h({{x}}_{li}| \lambda_{lk},\alpha_l).
#' }
#' The density functions (or probability distribution functions)
#' \deqn{h(.|\lambda_{lk},\alpha_l)}
#' can be any implemented model (Gaussian, Poisson,...).
#'
#' @slot ldata  a list of IClusterComponent.
#' @seealso [\code{\linkS4class{IClusterModelBase}}] class
#'
#' @examples
#' getSlots("ClusterHeterogeneous")
#'
#' @author Serge Iovleff
#'
#' @name ClusterHeterogeneous
#' @rdname ClusterHeterogeneous-class
#' @aliases ClusterHeterogeneous-class
#' @exportClass ClusterHeterogeneous
#'
setClass(
    Class="ClusterHeterogeneous",
    representation( ldata = "list"),
    contains=c("IClusterModelBase"),
    validity=function(object)
    {
      nbData = length(object@ldata)
      if (nbData == 0) {stop("At least on data set must be given.");}
      for (l in 1:nbData)
      {
        if (nrow(object@ldata[[1]]@data) != object@nbSample)
        {stop("All data sets must have the same number of individuals (number of rows).");}
      }
      return(TRUE)
    }
)

#' Initialize an instance of a MixAll S4 class.
#'
#' Initialization method of the [\code{\linkS4class{ClusterHeterogeneous}}] class.
#' Used internally in the 'MixAll' package.
#'
#' @rdname initialize-methods
#' @keywords internal
setMethod(
    f="initialize",
    signature=c("ClusterHeterogeneous"),
    definition=function(.Object, ldata =list(), nbCluster=2)
    {
      # for data
      if(missing(ldata)) {stop("ldata is mandatory in ClusterHeterogeneous.")}
      nbData = length(ldata)
      if (nbData == 0) {stop("At least on data set must be given.")}
      .Object@ldata <- ldata;
      # take first element of the list, this will give us the dimensions
      nbSample = nrow(.Object@ldata[[1]]@data);
      .Object <- callNextMethod(.Object, nbSample, nbCluster)
      # validate
      validObject(.Object)
      return(.Object)
    }
)

#' @rdname print-methods
#' @aliases print print,ClusterHeterogeneous-method
#'
setMethod(
  f="print",
  signature=c("ClusterHeterogeneous"),
  function(x,...){
    cat("****************************************\n")
    callNextMethod()
    nbData <- length(x@ldata)
    if(nbData>0)
    {
      for (l in 1:nbData)
      {
        cat("* model name = ", x@ldata[[l]]@modelName, "\n")
        print(format(x@ldata[[l]]@data),quote=FALSE);
      }
    }
    cat("****************************************\n")
    if(nbData>0)
    {
      for (l in 1:nbData)
      {
        for(k in 1:length(x@pk))
        {
          cat("*** Cluster: ",k,"\n")
          cat("* Proportion = ", format(x@pk[k]), "\n")
          cat("* model name = ", x@ldata[[l]]@modelName, "\n");
          print(x@ldata[[l]],k);
        }
      }
      cat("****************************************\n")
    }
  }
)

#' @rdname show-methods
#' @aliases show-ClusterHeterogeneous,ClusterHeterogeneous,ClusterHeterogeneous-method
setMethod(
    f="show",
    signature=c("ClusterHeterogeneous"),
    function(object)
    {
      cat("****************************************\n")
      callNextMethod()
      nbData <- length(object@ldata)
      if(nbData>0)
      {
        for (l in 1:nbData)
        {
          cat("* model name = ", object@ldata[[l]]@modelName, "\n")
          nrowShow <- min(10,nrow(object@ldata[[l]]@data))
          ncolShow <- min(10,ncol(object@ldata[[l]]@data))
          cat("* data (limited to 10 samples and 10 variables) =\n")
          print(format(object@ldata[[l]]@data[1:nrowShow,1:ncolShow]),quote=FALSE)
        }
      }
      cat("* ... ...\n")

      cat("****************************************\n")
      if(nbData>0)
      {
        for (l in 1:nbData)
        {
          for(k in 1:length(object@pk))
          {
            cat("*** Cluster: ",k,"\n")
            cat("* Proportion = ", format(object@pk[k]),"\n")
            cat("*\n")
            cat("* model name = ", object@ldata[[l]]@modelName, "\n");
            print(object@ldata[[l]], k);
          }
        }
        cat("****************************************\n")
      }
    }
)

#' @rdname summary-methods
#' @aliases summary summary,ClusterHeterogeneous-method
setMethod(
    f="summary",
    signature=c("ClusterHeterogeneous"),
    function(object, ...)
    {
      cat("**************************************************************\n")
      callNextMethod()
      cat("**************************************************************\n")
    }
)

#' Plotting of a class [\code{\linkS4class{ClusterHeterogeneous}}]
#'
#' Plotting data from a [\code{\linkS4class{ClusterHeterogeneous}}] object
#' using the estimated parameters and partition.
#'
#' @param x an object of class [\code{\linkS4class{ClusterHeterogeneous}}]
#' @param y a number between 1 and K-1.
#' @param ... further arguments passed to or from other methods
#'
#' @importFrom graphics plot
#' @aliases plot-ClusterHeterogeneous
#' @docType methods
#' @rdname plot-ClusterHeterogeneous-method
#' @export
#'
#' @seealso \code{\link{plot}}
#' @examples
#' \dontrun{
#'   ## the car data set
#'   data(car)
#'   model <- clusterHeterogeneous(car, 3, strategy = clusterFastStrategy())
#'   plot(model)
#'   }
#'
setMethod(
    f="plot",
    signature=c("ClusterHeterogeneous"),
    function(x, y, ...)
    {
      # total number of cluster in the data set
      nbCluster = ncol(x@tik);
      # check y, no y => display all dimensions
      if (missing(y)) { y=1:(nbCluster-1); }
      else
      { if (round(y)!=y) {stop("y must be an integer.")}
        if (y>nbCluster-1)
          stop("y should not be greater than K-1")
        y <- 1:y
      }
      # get representation
      Y=.visut(x@tik, nbCluster);
      if (nbCluster == 2) { ndim = 1;}
      else { ndim = ncol(Y);}
      # Compute gaussian statistics
      mean  <- matrix(0, nrow = x@nbCluster, ncol =ndim)
      sigma <- matrix(1, nrow = x@nbCluster, ncol =ndim)
      for (k in 1:nbCluster)
      {
        wcov = cov.wt(as.matrix(Y), x@tik[,k], method = "ML");
        mean[k,]  = wcov$center;
        sigma[k,] = sqrt(diag(wcov$cov))
      }
      # create gaussian model
      gauss<-new("ClusterDiagGaussian", Y, nbCluster = x@nbCluster)
      gauss@component@mean = mean
      gauss@component@sigma= sigma
      gauss@pk   = x@pk
      gauss@tik  = x@tik
      gauss@lnFi = x@lnFi
      gauss@zi   = x@zi
      #gauss@component@missing     = x@component@missing
      gauss@lnLikelihood = x@lnLikelihood
      gauss@criterion    = x@criterion
      gauss@nbFreeParameter = x@nbFreeParameter
      gauss@strategy        = x@strategy
      .clusterPlot(gauss, y, .dGauss,...);
    }
)

# get logistic representation
.visut <- function(t, gp)
{ m <- min(t[,gp]);
  if (m==0) t[,gp] = t[,gp] + 1e-30;
  return(scale(log(sweep(t,1,t[,gp],FUN="/")+ 1e-30), center=TRUE, scale=FALSE)[,-gp])
}

