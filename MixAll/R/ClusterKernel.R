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
#' Create an instance of the [\code{\linkS4class{ClusterKernel}}] class
#'
#' This function computes the optimal kernel mixture model according
#' to the [\code{criterion}] among the number of clusters given in
#' [\code{nbCluster}], using the strategy specified in [\code{strategy}].
#'
#' @param data frame or matrix containing the data. Rows correspond to observations
#' and columns correspond to variables.
#' @param dim integer giving the dimension of the Gaussian density. Default is 10.
#' @param nbCluster  [\code{\link{vector}}] listing the number of clusters to test.
#' @param modelNames [\code{\link{vector}}] of model names to run. By default only
#' "kernelGaussian_pk_s" is estimated. All the model names are given by the method
#' [\code{\link{clusterKernelNames}}].
#' @param kernelName string with a kernel name. Possible values:
#' "gaussian", "polynomial", "exponential". Default is "gaussian".
#' @param kernelParameters [\code{\link{vector}}] with the parameters of
#' the chosen kernel. Default is 1.
#' @param strategy a [\code{\linkS4class{ClusterStrategy}}] object containing
#' the strategy to run. [\code{\link{clusterStrategy}}]() method by default.
#' @param criterion character defining the criterion to select the best model.
#' The best model is the one with the lowest criterion value.
#' Possible values: "BIC", "AIC", "ICL". Default is "ICL".
#' @param nbCore integer defining the number of processor to use (default is 1, 0 for all).
#'
#' @note in the ClusterKernel instance returned by the function, data is the Gram matrix.
#' @examples
#' ## A quantitative example with the famous geyser data set
#' data(bullsEye)
#' ## estimate model (using fast strategy, results may be misleading)
#' model <- clusterKernel( data=bullsEye[,1:2], nbCluster=2:3
#'                       , modelNames= "kernelGaussian_pk_s"
#'                       , strategy = clusterFastStrategy()
#'                       )
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
#'
#' @return An instance of the [\code{\linkS4class{ClusterKernel}}] class.
#' @author Serge Iovleff
#' @export
#'
clusterKernel <- function( data, dim = 10, nbCluster=2
                         , modelNames= "kernelGaussian_pk_s"
                         , kernelName = "gaussian", kernelParameters = 1.
                         , strategy=clusterStrategy(), criterion="ICL"
                         , nbCore = 1)
{
  # check nbCluster
  nbClusterModel = length(nbCluster);
  nbClusterMin = min(nbCluster);
  nbClusterMax = max(nbCluster);
  if (nbClusterMin < 1)
  { stop("The number of clusters must be greater or equal to 1")}
  
  # check data
  if (missing(data))
  {stop("data is mandatory in clusterKernel")}
  data = as.matrix(data)
  if (nrow(data) <= 3*nbClusterMax)
  { stop("There is not enough individuals (rows) in the data set")}
  if (ncol(data) < 1) {stop("Error: empty data set")}
  # check dim
  if (dim < 1)
  { stop("The dimension must be greater or equal to 1")}

  # check criterion
  if(sum(criterion %in% c("BIC","AIC", "ICL")) != 1)
  { stop("criterion is not valid. See ?clusterKernel for the list of valid criterion")}

  # check modelNames
  if (!clusterValidKernelNames(modelNames))
  { stop("modelNames is not valid. See ?clusterKernelNames for the list of valid model names")}

  # check kernelName
  if(sum(kernelName %in% c("gaussian","polynomial", "exponential")) != 1)
  { stop("kernelName is not valid. See ?clusterKernel for the list of valid kernel name")}
  if (is.null(kernelParameters)) { kernelParameters = c(1)}

  # check strategy
  if(class(strategy)[1] != "ClusterStrategy")
  {stop("strategy is not a ClusterStrategy class (must be an instance of the class ClusterStrategy).")}
  validObject(strategy);

  # Create model
  model = new("ClusterKernel", data, dim, kernelName, kernelParameters)
  model@strategy = strategy;
  # start estimation of the models
  resFlag <- FALSE;
  if (length(nbCluster) >0)
  {
    # start estimation of the models
    resFlag = .Call("clusterKernelMixture", model, nbCluster, modelNames, strategy, criterion, nbCore, PACKAGE="MixAll")
  }
  if (resFlag != TRUE ) {cat("WARNING: An error occur during the clustering process");}
  model
}

#-----------------------------------------------------------------------
#' Definition of the [\code{\linkS4class{ClusterKernelComponent}}] class
#'
#' This class defines a kernel component of a mixture Model. It inherits
#' from [\code{\linkS4class{IClusterComponent}}].
#'
#' @slot dim  Vector with the dimension of the kth cluster
#' @slot sigma  Vector with the standard deviation in the kth cluster.
#'
#' @seealso [\code{\linkS4class{IClusterComponent}}] class
#'
#' @examples
#' getSlots("ClusterKernelComponent")
#'
#' @author Serge Iovleff
#'
#' @name ClusterKernelComponent
#' @rdname ClusterKernelComponent-class
#' @aliases ClusterKernelComponent-class
#' @exportClass ClusterKernelComponent
#'
setClass(
  Class="ClusterKernelComponent",
  representation( sigma = "vector", dim = "vector"),
  contains=c("IClusterComponent"),
  validity=function(object)
  {
    if ( length(object@sigma) == 0) { stop("sigma must be a vector.")}
    if ( length(object@dim) == 0)   { stop("dim must be a vector")}
    if (!clusterValidKernelNames(object@modelName))
    {stop("Invalid Kernel mixture model name.")}
    return(TRUE)
  }
)
#' Initialize an instance of a MixAll S4 class.
#'
#' Initialization method of the [\code{\linkS4class{ClusterKernelComponent}}] class.
#' Used internally in the 'MixAll' package.
#'
#' @rdname initialize-methods
#' @keywords internal
setMethod(
    f="initialize",
    signature=c("ClusterKernelComponent"),
    definition=function( .Object, data=matrix(nrow=0, ncol=0), dim =10, nbCluster=2, modelName="kernelGaussian_pk_sk")
    {
      # for data
      if(missing(data)) { {stop("There is not enough individuals (rows) in the data set")}}
      # check model name
      if (!clusterValidKernelNames(modelName)) { stop("modelName is invalid");}
      # call base class initialize
      .Object <- callNextMethod(.Object, data, modelName)
      # create slots
      nbVariable    = ncol(.Object@data);
      .Object@sigma = rep(1., nbCluster);
      .Object@dim   = rep(dim, nbCluster);
      # validate
      validObject(.Object)
      return(.Object)
    }
)

#' @rdname extract-methods
#' @aliases [,ClusterKernelComponent-method
setMethod(
    f="[",
    signature(x = "ClusterKernelComponent"),
    definition=function(x, i, j, drop)
    {
      if ( missing(j) )
      {
        switch(EXPR=i,
            "sigma"  = {return(x@sigma)},
            "dim"    = {return(x@dim)},
            stop("This attribute doesn't exist !")
        )
      }
      else
      {
        if (!is.numeric(j)) {stop("j must be an integer.")}
        if (round(j)!=j)    {stop("j must be an integer.")}
        switch(EXPR=i,
            "sigma"  = {return(x@sigma[j,])},
            "dim"    = {return(x@dim[j,])},
            stop("This attribute doesn't exist !")
        )
      }
    }
)

#' @rdname print-methods
#' @aliases print print,ClusterKernelComponent-method
#'
setMethod(
  signature=c("ClusterKernelComponent"),
  f="print",
  function(x,k,...)
  {
    cat("* sigma  = ", format(x@sigma[k]), "\n")
    cat("* dim    = ", format(x@dim[k]), "\n")
  }
)

#' @rdname show-methods
#' @aliases show-ClusterKernelComponent,ClusterKernelComponent,ClusterKernelComponent-method
setMethod(
  f="show",
  signature=c("ClusterKernelComponent"),
  function(object)
  {
    cat("* sigma  = ", format(object@sigma), "\n")
    cat("* dim    = ", format(object@dim), "\n")
  }
)

#-----------------------------------------------------------------------
#' Definition of the [\code{\linkS4class{ClusterKernel}}] class
#'
#' This class defines a Kernel mixture Model (KMM).
#'
#' This class inherits from the [\code{\linkS4class{IClusterModelBase}}] class.
#' A Kernel mixture  model is a mixture model of the form:
#' \deqn{
#'   f({x}|\boldsymbol{\theta})
#'   =\sum_{k=1}^K p_k \prod_{j=1}^d \phi(x_j;\sigma^2_{k})
#'    \quad x \in {R}^d.
#' }
#' Some constraints can be added to the variances in order to reduce the number
#' of parameters.
#'
#' @slot component  A [\code{\linkS4class{ClusterKernelComponent}}] with the
#' dim and standard deviation of the kernel mixture model.
#' @slot rawData  A matrix with the original data set
#' @slot kernelName string with the name of the kernel to use. Possible values:
#' "gaussian", "polynomial", "exponential". Default is "gaussian".
#' @slot kernelParameters vector with the parameters of the kernel.
#' @seealso [\code{\linkS4class{IClusterModelBase}}] class
#'
#' @examples
#' getSlots("ClusterKernel")
#' data(geyser)
#' new("ClusterKernel", data=geyser)
#'
#' @author Serge Iovleff
#'
#' @name ClusterKernel
#' @rdname ClusterKernel-class
#' @aliases ClusterKernel-class
#' @exportClass ClusterKernel
#'
setClass(
  Class="ClusterKernel",
  representation( component = "ClusterKernelComponent"
                , rawData = "matrix"
                , kernelName = "character"
                , kernelParameters = "vector"
                ),
  contains=c("IClusterModelBase"),
  validity=function(object)
  {
    if (length(object@component@dim)!=object@nbCluster)
    {stop("mean must have nbCluster length.")}
    if (length(object@component@sigma)!=object@nbCluster)
    {stop("sigma must have nbCluster length.")}
    if (!clusterValidKernelNames(object@component@modelName))
    {stop("Invalid kernel mixture model name.")}
    # check kernelName
    if(sum(object@kernelName %in% c("gaussian","polynomial", "exponential")) != 1)
    { stop("kernelName is not valid. See ?ClusterKernel for the list of valid kernel name")}
    return(TRUE)
  }
)

#' Initialize an instance of a MixAll S4 class.
#'
#' Initialization method of the [\code{\linkS4class{ClusterKernel}}] class.
#' Used internally in the 'MixAll' package.
#'
#' @rdname initialize-methods
#' @keywords internal
setMethod(
    f="initialize",
    signature=c("ClusterKernel"),
    definition=function(.Object, data=matrix(nrow=0, ncol=0), dim= 10
                       , kernelName = "gaussian", kernelParameters = NULL
                       , nbCluster=2, modelName="kernelGaussian_pk_sk")
    {
      # for data
      if(missing(data)) {stop("data is mandatory in ClusterKernel.")}
      .Object@rawData = as.matrix(data);
      # initialize fields
      .Object@kernelName = kernelName;
      if (is.null(kernelParameters)) {.Object@kernelParameters = 1;}
      else                           {.Object@kernelParameters = kernelParameters;}
      # initialize component
      .Object@component = new("ClusterKernelComponent", data, dim, nbCluster, modelName);
      .Object <- callNextMethod(.Object, nrow(.Object@component@data), nbCluster);
      # validate
      validObject(.Object);
      return(.Object);
    }
)

#' @rdname print-methods
#' @aliases print print,ClusterKernel-method
#'
setMethod(
  f="print",
  signature=c("ClusterKernel"),
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
#' @aliases show-ClusterKernel,ClusterKernel,ClusterKernel-method
setMethod(
  f="show",
  signature=c("ClusterKernel"),
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
#' @aliases summary summary,ClusterKernel-method
#'
setMethod(
  f="summary",
  signature=c("ClusterKernel"),
  function(object, ...)
  {
    cat("**************************************************************\n")
    callNextMethod()
    summary(object@component);
    cat("**************************************************************\n")
  }
)

#' Plotting of a class [\code{\linkS4class{ClusterKernel}}]
#'
#' Plotting data from a [\code{\linkS4class{ClusterKernel}}] object
#' using the estimated parameters and partition.
#'
#' @param x an object of class [\code{\linkS4class{ClusterKernel}}]
#' @param y a list of variables to plot (subset). Variables names or indices.
#' If missing all the variables are represented.
#' @param ... further arguments passed to or from other methods
#'
#' @importFrom graphics plot
#' @aliases plot-ClusterKernel
#' @docType methods
#' @rdname plot-ClusterKernel-method
#' @export
#'
#' @seealso \code{\link{plot}}
#' @examples
#'  ## the famous iris data set
#' \dontrun{
#'   data(bullsEye)
#'   model <- clusterKernel( bullsEye[,1:2]
#'                         , 2
#'                         , modelNames= "kernelGaussian_pk_s"
#'                         , strategy = clusterFastStrategy()
#'                         )
#'   plot(model)
#'   }
#'
setMethod(
    f="plot",
    signature=c("ClusterKernel"),
    function(x, y, ...)
    { # use generic plot
      # total number of variable in the data set
      nbVariable = ncol(x@rawData);
      # no y => display all variables
      if (missing(y)) { y=1:nbVariable; }
      else # perform some check
      {
        if (is.numeric(y)) # numbers of the columns to plot are given
        {
          if (max(y)>nbVariable)
            stop("In .clusterPlot, y indices mismatch the data dimension")
        }
        else # names of the variables to plot are given
        {
          if ( sum(y %in% colnames(x@rawData))!= length(y) )
          { stop(cat("In plot, unknown variables: ", paste(y[which(!(y %in% colnames(x@rawData)))]),"\n"))}
        }
      }
      # scatter plot
      plot(x@rawData[,y], col = x@zi+2)
    }
)

