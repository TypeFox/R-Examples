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
#' Create an instance of the [\code{\linkS4class{ClusterGamma}}] class
#'
#' This function computes the optimal gamma mixture model according
#' to the \code{criterion} among the list of model given in \code{modelNames}
#' and the number of clusters given in \code{nbCluster}, using the strategy
#' specified in \code{strategy}.
#'
#' @param data frame or matrix containing the data. Rows correspond to observations
#' and columns correspond to variables. If the data set contains NA values, they
#' will be estimated during the estimation process.
#' @param nbCluster  [\code{\link{vector}}] listing the number of clusters to test.
#' @param modelNames [\code{\link{vector}}] of model names to run. By default all
#' gamma models with free shape are estimated.  All the model names are given by
#' the method [\code{\link{clusterGammaNames}}].
#' @param strategy a [\code{\linkS4class{ClusterStrategy}}] object containing
#' the strategy to run. [\code{\link{clusterStrategy}}]() method by default.
#' @param criterion character defining the criterion to select the best model.
#' The best model is the one with the lowest criterion value.
#' Possible values: "BIC", "AIC", "ICL". Default is "ICL".
#' @param nbCore integer defining the number of processor to use (default is 1, 0 for all).
#'
#' @examples
#' ## A quantitative example with the famous geyser data set
#' data(geyser)
#' ## add 10 missing values
#' x = geyser;
#' x[round(runif(5,1,nrow(geyser))), 1] <- NA
#' x[round(runif(5,1,nrow(geyser))), 2] <- NA
#'
#' ## use graphics functions
#' model <- clusterGamma( data=x, nbCluster=2:3
#'                      , modelNames="gamma_pk_ajk_bjk"
#'                      , strategy = clusterFastStrategy())
#'
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
#' @return An instance of the [\code{\linkS4class{ClusterGamma}}] class.
#' @author Serge Iovleff
#' @export
#'
clusterGamma <- function( data, nbCluster=2
                        , modelNames= "gamma_pk_ajk_bjk"
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
  if(sum(criterion %in% c("BIC","AIC","ICL")) != 1)
  { stop("criterion is not valid. See ?clusterGamma for the list of valid criterion")}

  # check data
  data = as.matrix(data)
  if (nrow(data) <= 3*nbClusterMax) {stop("There is not enough individuals (rows) in the data set")}
  if (ncol(data) < 1) {stop("Error: empty data set")}

  # check modelNames
  if (!clusterValidGammaNames(modelNames))
  { stop("modelNames is not valid. See ?clusterGammaNames for the list of valid model names")}

  # check strategy
  if(class(strategy)[1] != "ClusterStrategy")
  {stop("strategy is not a Cluster Stategy class (must be an instance of the class ClusterStrategy).")}
  validObject(strategy);

  # Create model
  model = new("ClusterGamma", data)
  model@strategy = strategy;

  # start estimation of the models
  resFlag = .Call("clusterMixture", model, nbCluster, modelNames, strategy, criterion, nbCore, PACKAGE="MixAll")

  # set names
  if (resFlag != 1) {cat("WARNING: An error occur during the clustering process")}
  colnames(model@component@shape) <- colnames(model@component@data)
  colnames(model@component@scale) <- colnames(model@component@data)
  model
}

#-----------------------------------------------------------------------
#' Definition of the [\code{\linkS4class{ClusterGammaComponent}}] class
#'
#' This class defines a gamma component of a mixture Model. It inherits
#' from [\code{\linkS4class{IClusterComponent}}].
#'
#' @slot shape  Matrix with the shapes of the jth variable in the kth cluster.
#' @slot scale  Matrix with the scales of the jth variable in the kth cluster.
#'
#' @seealso [\code{\linkS4class{IClusterComponent}}] class
#'
#' @examples
#' getSlots("ClusterGammaComponent")
#'
#' @author Serge Iovleff
#'
#' @name ClusterGammaComponent
#' @rdname ClusterGammaComponent-class
#' @aliases ClusterGammaComponent-class
#' @exportClass ClusterGammaComponent
#'
setClass(
    Class="ClusterGammaComponent",
    representation( shape = "matrix", scale = "matrix"),
    contains=c("IClusterComponent"),
    validity=function(object)
    {
      if (ncol(object@shape)!=ncol(object@data))
      {stop("shape must have nbVariable columns.")}
      if (ncol(object@scale)!=ncol(object@data))
      {stop("scale must have nbVariable columns.")}
      if (!clusterValidGammaNames(object@modelName))
      {stop("Invalid Gaussian model name.")}
      return(TRUE)
    }
)
#' Initialize an instance of a MixAll S4 class.
#'
#' Initialization method of the [\code{\linkS4class{ClusterGammaComponent}}] class.
#' Used internally in the 'MixAll' package.
#'
#' @rdname initialize-methods
#' @keywords internal
setMethod(
    f="initialize",
    signature=c("ClusterGammaComponent"),
    definition=function(.Object, data=matrix(nrow=0, ncol=0), nbCluster=2, modelName="gamma_pk_ajk_bjk")
    {
      # for data
      if(missing(data)) {stop("data is mandatory in ClusterGammaComponent.")}
      # check model name
      if (!clusterValidGammaNames(modelName)) { stop("modelName is invalid");}
      # call base class
      .Object <- callNextMethod(.Object, data, modelName)
      # create slots
      nbVariable = ncol(.Object@data);
      .Object@shape = matrix(1., nrow=nbCluster, ncol=nbVariable);
      .Object@scale = matrix(1., nrow=nbCluster, ncol=nbVariable);
      # validate
      validObject(.Object)
      return(.Object)
    }
)

#' @rdname extract-methods
#' @aliases [,ClusterGammaComponent-method
setMethod(
  f="[",
  signature(x = "ClusterGammaComponent"),
  definition=function(x, i, j, drop)
  {
    if ( missing(j) )
    {
      switch(EXPR=i,
             "shape"={return(x@shape)},
             "scale"={return(x@scale)},
             stop("This attribute doesn't exist !")
            )
    }
    else
    {
      if (!is.numeric(j)) {stop("j must be an integer.")}
      if (round(j)!=j)    {stop("j must be an integer.")}
      switch(EXPR=i,
             "shape"={return(x@shape[j,])},
             "scale"={return(x@scale[j,])},
             stop("This attribute doesn't exist !")
            )
    }
  }
)

#' @rdname print-methods
#' @aliases print print,ClusterGammaComponent-method
#'
setMethod(
  f="print",
  signature=c("ClusterGammaComponent"),
  function(x,k,...)
  {
    cat("* Shapes     = ", format(x@shape[k,]), "\n")
    cat("* Scales     = ", format(x@scale[k,]), "\n")
  }
)

#' @rdname show-methods
#' @aliases show-ClusterGammaComponent,ClusterGammaComponent,ClusterGammaComponent-method
setMethod(
  f="show",
  signature=c("ClusterGammaComponent"),
  function(object)
  {
    cat("* Shapes     = ", format(object@shape), "\n")
    cat("* Scales     = ", format(object@scale), "\n")
  }
)

#-----------------------------------------------------------------------
#' Definition of the [\code{\linkS4class{ClusterGamma}}] class
#'
#' This class inherits from the [\code{\linkS4class{IClusterModelBase}}] class.
#' A gamma mixture model is a mixture model of the form:
#' \deqn{
#'   f({x}|\boldsymbol{\theta}) \\
#'   =\sum_{k=1}^K p_k \prod_{j=1}^d \gamma(x_j;a_{jk},b_{jk}) \\
#'    \quad {x} \in {R}^d.
#' }
#' Constraints can be added to the shapes and/or scales in order to reduce the number
#' of parameters.
#'
#' @slot component  A [\code{\linkS4class{ClusterGammaComponent}}] with the
#' shapes and the scales of the component mixture model.
#' @seealso [\code{\linkS4class{IClusterModelBase}}] class
#'
#' @examples
#'   getSlots("ClusterGamma")
#'   data(geyser)
#'   new("ClusterGamma", data=geyser)
#'
#' @author Serge Iovleff
#'
#' @name ClusterGamma
#' @rdname ClusterGamma-class
#' @aliases ClusterGamma-class
#' @exportClass ClusterGamma
#'
setClass(
    Class="ClusterGamma",
    representation( component = "ClusterGammaComponent"),
    contains=c("IClusterModelBase"),
    validity=function(object)
    {
      if (ncol(object@component@shape)!=ncol(object@component@data))
      {stop("shape must have nbVariable columns.")}
      if (nrow(object@component@shape)!=object@nbCluster)
      {stop("shape must have nbCluster rows.")}

      if (ncol(object@component@scale)!=ncol(object@component@data))
      {stop("scale must have nbVariable columns.")}
      if (nrow(object@component@scale)!=object@nbCluster)
      {stop("scale must have nbCluster rows.")}

      if (!clusterValidGammaNames(object@component@modelName))
      {stop("Invalid Gamma model name.")}
      return(TRUE)
    }
)

#' Initialize an instance of a MixAll S4 class.
#'
#' Initialization method of the [\code{\linkS4class{ClusterGamma}}] class.
#' Used internally in the 'MixAll' package.
#'
#' @rdname initialize-methods
#' @keywords internal
setMethod(
    f="initialize",
    signature=c("ClusterGamma"),
    definition=function(.Object, data=matrix(nrow=0, ncol=0), nbCluster=2, modelName="gamma_pk_ajk_bjk")
    {
      # for data
      if(missing(data)) {stop("data is mandatory in ClusterGamma.")}
      .Object@component = new("ClusterGammaComponent", data, nbCluster, modelName);
      .Object <- callNextMethod(.Object, nrow(.Object@component@data), nbCluster);
      # validate
      validObject(.Object)
      return(.Object)
    }
)

#' @rdname print-methods
#' @aliases print print,ClusterGamma-method
#'
setMethod(
  f="print",
  signature=c("ClusterGamma"),
  function(x,...){
    cat("****************************************\n")
    callNextMethod();
    cat("****************************************\n")
    for(k in 1:length(x@pk))
    {
      cat("*** Cluster: ",k,"\n");
      cat("* Proportion = ", format(x@pk[k]), "\n");
      print(x@component, k);
      cat("****************************************\n")
    }
  }
)

#' @rdname show-methods
#' @aliases show-ClusterGamma,ClusterGamma,ClusterGamma-method
setMethod(
  f="show",
  signature=c("ClusterGamma"),
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
#' @aliases summary summary,MixmodResults-method
#'
setMethod(
  f="summary",
  signature=c("ClusterGamma"),
  function(object, ...)
  {
    cat("**************************************************************\n")
    callNextMethod()
    summary(object@component);
    cat("**************************************************************\n")
  }
)

#' Plotting of a class [\code{\linkS4class{ClusterGamma}}]
#'
#' Plotting data from a [\code{\linkS4class{ClusterGamma}}] object
#' using the estimated parameters and partition.
#'
#' @param x an object of class [\code{\linkS4class{ClusterGamma}}]
#' @param y a list of variables to plot (subset). Variables names or indices.
#' If missingValues all the variables are represented.
#' @param ... further arguments passed to or from other methods
#'
#' @importFrom graphics plot
#' @aliases plot-ClusterGamma, ClusterGamma-method
#' @docType methods
#' @rdname plot-ClusterGamma-method
#' @export
#'
#' @seealso \code{\link{plot}}
#' @examples
#'   ## for quantitative case
#' \dontrun{
#'   data(iris)
#'   model <- clusterGamma( data=iris[1:4], nbCluster=3
#'                        , modelNames=clusterGammaNames(prop = "equal")
#'                        , strategy = clusterFastStrategy())
#'   plot(model)
#'   plot(model, c(1,3))
#'   plot(model, c("Sepal.Length","Sepal.Width"))
#'   }
#'
setMethod(
    f="plot",
    signature=c("ClusterGamma"),
    function(x, y, ...)
    { # use generic plot
      .clusterPlot(x, y, .dGamma,...);
    }
)

# wrapper of dnorm
# x a vector with the point
.dGamma <- function(x, j, k, model)
{ dgamma(x, shape = (model@component@shape)[k, j] , scale = (model@component@scale)[k, j])}
