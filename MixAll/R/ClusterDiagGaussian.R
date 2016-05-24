#-----------------------------------------------------------------------
#     Copyright (C) 2012-2015  Serge Iovleff, University Lille 1, Inria
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
#' Create an instance of the [\code{\linkS4class{ClusterDiagGaussian}}] class
#'
#' This function computes the optimal diagonal Gaussian mixture model according
#' to the \code{criterion} among the list of model given in \code{modelNames}
#' and the number of clusters given in \code{nbCluster}, using the strategy
#' specified in \code{strategy}.
#'
#' @param data frame or matrix containing the data. Rows correspond to observations
#' and columns correspond to variables. If the data set contains NA values, they
#' will be estimated during the estimation process.
#' @param nbCluster  [\code{\link{vector}}] listing the number of clusters to test.
#' @param modelNames [\code{\link{vector}}] of model names to run. By default all diagonal
#' Gaussian models are estimated. All the model names are given by the method
#' [\code{\link{clusterDiagGaussianNames}}].
#' @param strategy a [\code{\linkS4class{ClusterStrategy}}] object containing
#' the strategy to run. [\code{\link{clusterStrategy}}]() method by default.
#' @param criterion character defining the criterion to select the best model.
#' The best model is the one with the lowest criterion value.
#' Possible values: "BIC", "AIC", "ICL". Default is "ICL".
#' @param nbCore integer defining the number of processors to use (default is 1, 0 for all).
#'
#' @examples
#' ## A quantitative example with the famous geyser data set
#' data(geyser)
#' ## add 10 missing values as random
#' x = as.matrix(geyser); n <- nrow(x); p <- ncol(x);
#' indexes <- matrix(c(round(runif(5,1,n)), round(runif(5,1,p))), ncol=2);
#' x[indexes] <- NA;
#' ## estimate model (using fast strategy, results may be misleading)
#' model <- clusterDiagGaussian( data=x, nbCluster=2:3
#'                             , modelNames=c( "gaussian_pk_sjk")
#'                             , strategy = clusterFastStrategy()
#'                             )
#'
#' ## use graphics functions
#' \dontrun{
#' plot(model)
#' }
#'
#' ## get summary
#' summary(model)
#' ## print model
#' \dontrun{
#' print(model)
#' }
#' ## get estimated missing values
#' missingValues(model)
#'
#' @return An instance of the [\code{\linkS4class{ClusterDiagGaussian}}] class.
#' @author Serge Iovleff
#' @export
#'
clusterDiagGaussian <- function( data, nbCluster=2
                               , modelNames=clusterDiagGaussianNames()
                               , strategy=clusterFastStrategy()
                               , criterion="ICL"
                               , nbCore = 1)
{
  # check nbCluster
  nbClusterModel = length(nbCluster);
  nbClusterMin = min(nbCluster);
  nbClusterMax = max(nbCluster);
  if (nbClusterMin < 1) { stop("The number of clusters must be greater or equal to 1")}

  # check criterion
  if(sum(criterion %in% c("BIC","AIC", "ICL")) != 1)
  { stop("criterion is not valid. See ?clusterDiagGaussian for the list of valid criterion")}

  # check data
  data = as.matrix(data)
  if (nrow(data) <= 3*nbClusterMax) {stop("There is not enough individuals (rows) in the data set")}
  if (ncol(data) < 1) {stop("Error: empty data set")}

  # check modelNames
  if (!clusterValidDiagGaussianNames(modelNames))
  { stop("modelNames is not valid. See ?clusterDiagGaussianNames for the list of valid model names")}

  # check strategy
  if(class(strategy)[1] != "ClusterStrategy")
  {stop("strategy is not a Cluster Strategy class (must be an instance of the class ClusterStrategy).")}
  validObject(strategy);

  # Create model
  model = new("ClusterDiagGaussian", data)
  model@strategy = strategy;
  # start estimation of the models
  resFlag = .Call("clusterMixture", model, nbCluster, modelNames, strategy, criterion, nbCore, PACKAGE="MixAll");
  # set names
  if (resFlag != TRUE ) {cat("WARNING: An error occur during the clustering process");}
  colnames(model@component@mean)  <- colnames(model@component@data);
  colnames(model@component@sigma) <- colnames(model@component@data);
  model
}

#-----------------------------------------------------------------------
#' Definition of the [\code{\linkS4class{ClusterDiagGaussianComponent}}] class
#'
#' This class defines a diagonal Gaussian component of a mixture model. It inherits
#' from [\code{\linkS4class{IClusterComponent}}].
#'
#' @slot mean  Matrix with the mean of the jth variable in the kth cluster.
#' @slot sigma  Matrix with the standard deviation of the jth variable in the kth cluster.
#'
#' @seealso [\code{\linkS4class{IClusterComponent}}] class
#'
#' @examples
#' getSlots("ClusterDiagGaussianComponent")
#'
#' @author Serge Iovleff
#'
#' @name ClusterDiagGaussianComponent
#' @rdname ClusterDiagGaussianComponent-class
#' @aliases ClusterDiagGaussianComponent-class
#' @exportClass ClusterDiagGaussianComponent
#'
setClass(
  Class="ClusterDiagGaussianComponent",
  representation( mean = "matrix", sigma = "matrix"),
  contains=c("IClusterComponent"),
  validity=function(object)
  {
    if (ncol(object@mean)!=ncol(object@data))
    {stop("mean must have nbVariable columns.")}
    if (ncol(object@sigma)!=ncol(object@data))
    {stop("sigma must have nbVariable columns.")}
    if (!clusterValidDiagGaussianNames(object@modelName))
    {stop("Invalid Gaussian model name.")}
    return(TRUE)
  }
)
#' Initialize an instance of a MixAll S4 class.
#'
#' Initialization method of the [\code{\linkS4class{ClusterDiagGaussianComponent}}] class.
#' Used internally in the 'MixAll' package.
#'
#' @rdname initialize-methods
#' @keywords internal
setMethod(
    f="initialize",
    signature=c("ClusterDiagGaussianComponent"),
    definition=function(.Object, data=matrix(nrow=0, ncol=0), nbCluster=2, modelName="gaussian_pk_sjk")
    {
      # for data
      if(missing(data)) {stop("data is mandatory in ClusterDiagGaussianComponent.")}
      # check model name
      if (!clusterValidDiagGaussianNames(modelName)) { stop("modelName is invalid");}
      # call base class initialize
      .Object <- callNextMethod(.Object, data, modelName)
      # create slots
      nbVariable = ncol(.Object@data);
      .Object@mean  = matrix(0., nrow=nbCluster, ncol=nbVariable);
      .Object@sigma = matrix(1., nrow=nbCluster, ncol=nbVariable);
      # validate
      validObject(.Object)
      return(.Object)
    }
)

#' @rdname extract-methods
#' @aliases [,ClusterDiagGaussianComponent-method
setMethod(
    f="[",
    signature(x = "ClusterDiagGaussianComponent"),
    definition=function(x, i, j, drop)
    {
      if ( missing(j) )
      {
        switch(EXPR=i,
            "mean" ={return(x@mean)},
            "sigma"={return(x@sigma)},
            stop("This attribute doesn't exist !")
        )
      }
      else
      {
        if (!is.numeric(j)) {stop("j must be an integer.")}
        if (round(j)!=j)    {stop("j must be an integer.")}
        switch(EXPR=i,
            "mean" ={return(x@mean[j,])},
            "sigma"={return(x@sigma[j,])},
            stop("This attribute doesn't exist !")
        )
      }
    }
)

#' @rdname print-methods
#' @aliases print print,ClusterDiagGaussianComponent-method
#'
setMethod(
  signature=c("ClusterDiagGaussianComponent"),
  f="print",
  function(x,k,...)
  {
    cat("* Means     = ", format(x@mean[k,]), "\n")
    cat("* S.D.      = ", format(x@sigma[k,]), "\n")
  }
)

#' @rdname show-methods
#' @aliases show-ClusterDiagGaussianComponent,ClusterDiagGaussianComponent,ClusterDiagGaussianComponent-method
setMethod(
  f="show",
  signature=c("ClusterDiagGaussianComponent"),
  function(object)
  {
    cat("* Means     = ", format(object@mean), "\n")
    cat("* S.D.      = ", format(object@sigma), "\n")
  }
)


#-----------------------------------------------------------------------
#' Definition of the [\code{\linkS4class{ClusterDiagGaussian}}] class
#'
#' This class defines a diagonal Gaussian mixture Model.
#'
#' This class inherits from the [\code{\linkS4class{IClusterModelBase}}] class.
#' A diagonal gaussian model is a mixture model of the form:
#' \deqn{
#'   f({x}|\boldsymbol{\theta})
#'   =\sum_{k=1}^K p_k \prod_{j=1}^d \phi(x_j;\mu_{jk},\sigma^2_{jk})
#'    \quad x \in {R}^d.
#' }
#' Some constraints can be added to the variances in order to reduce the number
#' of parameters.
#'
#' @slot component  A [\code{\linkS4class{ClusterDiagGaussianComponent}}] with the
#' mean and standard deviation of the diagonal mixture model.
#' @seealso [\code{\linkS4class{IClusterModelBase}}] class
#'
#' @examples
#' getSlots("ClusterDiagGaussian")
#' data(geyser)
#' new("ClusterDiagGaussian", data=geyser)
#'
#' @author Serge Iovleff
#'
#' @name ClusterDiagGaussian
#' @rdname ClusterDiagGaussian-class
#' @aliases ClusterDiagGaussian-class
#' @exportClass ClusterDiagGaussian
#'
setClass(
  Class="ClusterDiagGaussian",
  representation( component = "ClusterDiagGaussianComponent"),
  contains=c("IClusterModelBase"),
  validity=function(object)
  {
    if (nrow(object@component@mean)!=object@nbCluster)
    {stop("mean must have nbCluster rows.")}
    if (ncol(object@component@mean)!=ncol(object@component@data))
    {stop("mean must have nbVariable columns.")}
    if (nrow(object@component@sigma)!=object@nbCluster)
    {stop("sigma must have nbCluster rows.")}
    if (ncol(object@component@sigma)!=ncol(object@component@data))
    {stop("sigma must have nbVariable columns.")}
    if (!clusterValidDiagGaussianNames(object@component@modelName))
    {stop("Invalid Gaussian model name.")}
    return(TRUE)
  }
)

#' Initialize an instance of a MixAll S4 class.
#'
#' Initialization method of the [\code{\linkS4class{ClusterDiagGaussian}}] class.
#' Used internally in the 'MixAll' package.
#'
#' @rdname initialize-methods
#' @keywords internal
setMethod(
    f="initialize",
    signature=c("ClusterDiagGaussian"),
    definition=function(.Object, data=matrix(nrow=0, ncol=0), nbCluster=2, modelName="gaussian_pk_sjk")
    {
      # for data
      if(missing(data)) {stop("data is mandatory in ClusterDiagGaussian.")}
      # initialize component
      .Object@component = new("ClusterDiagGaussianComponent", data, nbCluster, modelName);
      .Object <- callNextMethod(.Object, nrow(.Object@component@data), nbCluster);
      # validate
      validObject(.Object)
      return(.Object)
    }
)

#' @rdname print-methods
#' @aliases print print,ClusterDiagGaussian-method
#'
setMethod(
  f="print",
  signature=c("ClusterDiagGaussian"),
  function(x,...){
    cat("****************************************\n")
    callNextMethod();
    cat("****************************************\n")
    for(k in 1:length(x@pk))
    {
      cat("*** Cluster: ",k,"\n")
      cat("* Proportion = ", format(x@pk[k]), "\n")
      print(x@component, k);
      cat("****************************************\n")
    }
  }
)

#' @rdname show-methods
#' @aliases show-ClusterDiagGaussian,ClusterDiagGaussian,ClusterDiagGaussian-method
setMethod(
  f="show",
  signature=c("ClusterDiagGaussian"),
  function(object)
  {
    cat("****************************************\n")
    callNextMethod();
    show(object@component);
    cat("****************************************\n")
    for(k in 1:length(object@pk))
    {
      cat("*** Cluster: ",k,"\n")
      cat("* Proportion = ", format(object@pk[k]), "\n")
      print(object@component, k);
      cat("****************************************\n")
    }
  }
)

#' @rdname summary-methods
#' @aliases summary summary,ClusterDiagGaussian-method
#'
setMethod(
  f="summary",
  signature=c("ClusterDiagGaussian"),
  function(object, ...)
  {
    cat("**************************************************************\n")
    callNextMethod()
    summary(object@component);
    cat("**************************************************************\n")
  }
)

#' Plotting of a class [\code{\linkS4class{ClusterDiagGaussian}}]
#'
#' Plotting data from a [\code{\linkS4class{ClusterDiagGaussian}}] object
#' using the estimated parameters and partition.
#'
#' @param x an object of class [\code{\linkS4class{ClusterDiagGaussian}}]
#' @param y a list of variables to plot (subset). Variables names or indices.
#' If missing all the variables are represented.
#' @param ... further arguments passed to or from other methods
#'
#' @importFrom graphics plot
#' @aliases plot-ClusterDiagGaussian
#' @docType methods
#' @rdname plot-ClusterDiagGaussian-method
#' @export
#'
#' @seealso \code{\link{plot}}
#' @examples
#'   ## the famous iris data set
#' \dontrun{
#'   data(iris)
#'   model <- clusterDiagGaussian(iris[1:4], 3, strategy = clusterFastStrategy())
#'   plot(model)
#'   plot(model, c(1,3))
#'   plot(model, c("Sepal.Length","Sepal.Width"))
#'   }
#'
setMethod(
    f="plot",
    signature=c("ClusterDiagGaussian"),
    function(x, y, ...)
    { # use generic plot
      .clusterPlot(x, y, .dGauss,...);
    }
)

# wrapper of dnorm
# x a vector with the point
.dGauss <- function(x, j, k, model)
{ dnorm(x, (model@component@mean)[k, j] , (model@component@sigma)[k, j])}
