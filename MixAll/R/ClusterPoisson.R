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
#' Create an instance of the [\code{\linkS4class{ClusterPoisson}}] class
#'
#' This function computes the optimal poisson mixture model according
#' to the [\code{criterion}] among the list of model given in [\code{modelNames}]
#' and the number of clusters given in [\code{nbCluster}], using the strategy
#' specified in [\code{strategy}].
#'
#' @param data frame or matrix containing the data. Rows correspond to observations
#' and columns correspond to variables. If the data set contains NA values, they
#' will be estimated during the estimation process.
#' @param nbCluster  [\code{\link{vector}}] listing the number of clusters to test.
#' @param modelNames [\code{\link{vector}}] of model names to run. By default all
#' poisson models are estimated.  All the model names are given by
#' the method [\code{\link{clusterPoissonNames}}].
#' @param strategy a [\code{\linkS4class{ClusterStrategy}}] object containing
#' the strategy to run. [\code{\link{clusterStrategy}}]() method by default.
#' @param criterion character defining the criterion to select the best model.
#' The best model is the one with the lowest criterion value.
#' Possible values: "BIC", "AIC", "ICL". Default is "ICL".
#' @param nbCore integer defining the number of processor to use (default is 1, 0 for all).
#'
#' @examples
#' ## A quantitative example with the DebTrivedi data set.
#' data(DebTrivedi)
#' dt <- DebTrivedi[1:500, c(1, 6,8, 15)]
#'
#' model <- clusterPoisson( data=dt, nbCluster=2
#'                        , modelNames=clusterPoissonNames(prop = "equal")
#'                        , strategy = clusterFastStrategy())
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
#' @return An instance of the [\code{\linkS4class{ClusterPoisson}}] class.
#' @author Serge Iovleff
#' @export
#'
clusterPoisson <- function( data, nbCluster=2
                          , modelNames= clusterPoissonNames()
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
  { stop("criterion is not valid. See ?clusterPoisson for the list of valid criterion")}

  # check data
  data = as.matrix(data)
  if (nrow(data) <= 3*nbClusterMax) {stop("There is not enough individuals (rows) in the data set")}
  if (ncol(data) < 1) {stop("Error: empty data set")}

  # check modelNames
  if (is.null(modelNames)) { modelNames = clusterPoissonNames()}
  if (!clusterValidPoissonNames(modelNames))
  { stop("modelNames is not valid. See ?clusterPoissonNames for the list of valid model names")}

  # check strategy
  if(class(strategy)[1] != "ClusterStrategy")
  {stop("strategy is not a Cluster Stategy class (must be an instance of the class ClusterStrategy).")}
  validObject(strategy);

  # Create model
  model = new("ClusterPoisson", data)
  model@strategy = strategy;

  # start estimation of the models
  resFlag = .Call("clusterMixture", model, nbCluster, modelNames, strategy, criterion, nbCore, PACKAGE="MixAll")

  # set names
  if (resFlag != 1) {cat("WARNING: An error occur during the clustering process")}
  colnames(model@component@lambda) <- colnames(model@component@data)
  model
}

#-----------------------------------------------------------------------
#' Definition of the [\code{\linkS4class{ClusterPoissonComponent}}] class
#'
#' This class defines a poisson component of a mixture Model. It inherits
#' from [\code{\linkS4class{IClusterComponent}}].
#'
#' @slot lambda Matrix with the mean of the jth variable in the kth cluster.
#' @seealso [\code{\linkS4class{IClusterComponent}}] class
#'
#' @examples
#' getSlots("ClusterPoissonComponent")
#'
#' @author Serge Iovleff
#'
#' @name ClusterPoissonComponent
#' @rdname IClusterComponent-class
#' @aliases ClusterPoissonComponent-class
#' @exportClass ClusterPoissonComponent
#'
setClass(
    Class="ClusterPoissonComponent",
    representation( lambda = "matrix"),
    contains=c("IClusterComponent"),
    validity=function(object)
    {
      if (ncol(object@lambda)!=ncol(object@data))
      {stop("lambda must have nbVariable columns.")}
      if (!clusterValidPoissonNames(object@modelName))
      {stop("Invalid poisson model name.")}
      return(TRUE)
    }
)
#' Initialize an instance of a MixAll S4 class.
#'
#' Initialization method of the [\code{\linkS4class{ClusterPoissonComponent}}] class.
#' Used internally in the 'MixAll' package.
#'
#' @rdname initialize-methods
#' @keywords internal
setMethod(
    f="initialize",
    signature=c("ClusterPoissonComponent"),
    definition=function(.Object, data=matrix(nrow=0, ncol=0), nbCluster=2, modelName="poisson_pk_ljk")
    {
      # for data
      if(missing(data)) {stop("data is mandatory in ClusterPoissonComponent.")}
      # check model name
      if (!clusterValidPoissonNames(modelName))
      { stop(paste(modelName,"modelName is invalid in ClusterPoissonComponent."));}
      # call base class
      .Object <- callNextMethod(.Object, data, modelName)
      # create slots
      nbVariable = ncol(.Object@data);
      .Object@lambda = matrix(1., nrow=nbCluster, ncol=nbVariable);
      # validate
      validObject(.Object)
      return(.Object)
    }
)

#' @rdname extract-methods
#' @aliases [,ClusterPoissonComponent-method
setMethod(
  f="[",
  signature(x = "ClusterPoissonComponent"),
  definition=function(x, i, j, drop)
  {
    if ( missing(j) )
    {
      switch(EXPR=i,
             "lambda"={return(x@lambda)},
             stop("This attribute doesn't exist !")
            )
    }
    else
    {
      if (!is.numeric(j)) {stop("j must be an integer.")}
      if (round(j)!=j)    {stop("j must be an integer.")}
      switch(EXPR=i,
             "lambda"={return(x@lambda[j,])},
             stop("This attribute doesn't exist !")
            )
    }
  }
)

#' @rdname print-methods
#' @aliases print print,ClusterPoissonComponent-method
#'
setMethod(
  f="print",
  signature=c("ClusterPoissonComponent"),
  function(x,k,...)
  {
    cat("* lambda     = ", format(x@lambda[k,]), "\n")
  }
)

#' @rdname show-methods
#' @aliases show-ClusterPoissonComponent,ClusterPoissonComponent,ClusterPoissonComponent-method
setMethod(
  f="show",
  signature=c("ClusterPoissonComponent"),
  function(object)
  {
    cat("* lambda     = ", format(object@lambda), "\n")
  }
)

#' Definition of the [\code{\linkS4class{ClusterPoisson}}] class
#'
#' This class inherits from the [\code{\linkS4class{IClusterModelBase}}] class.
#' A poisson mixture model is a mixture model of the form:
#' \deqn{
#'   f({x}|\boldsymbol{\theta}) \\
#'   =\sum_{k=1}^K p_k \prod_{j=1}^d Pois(x_j;\lambda_{jk}) \\
#'    \quad {x} \in {N}^d.
#' }
#'
#' @slot component  A [\code{\linkS4class{ClusterPoissonComponent}}] with the
#' lambda of the component mixture model.
#' @seealso [\code{\linkS4class{IClusterModelBase}}] class
#'
#' @examples
#'   getSlots("ClusterPoisson")
#'   data(DebTrivedi)
#'   dt <- DebTrivedi[, c(1, 6,8, 15)]
#'   new("ClusterPoisson", data=dt)
#'
#' @author Serge Iovleff
#'
#' @name ClusterPoisson
#' @rdname ClusterPoisson-class
#' @aliases ClusterPoisson-class
#' @exportClass ClusterPoisson
#'
setClass(
    Class="ClusterPoisson",
    representation( component = "ClusterPoissonComponent"),
    contains=c("IClusterModelBase"),
    validity=function(object)
    {
      if (ncol(object@component@lambda)!=ncol(object@component@data))
      {stop("lambda must have nbVariable columns.")}
      if (nrow(object@component@lambda)!=object@nbCluster)
      {stop("lambda must have nbCluster rows.")}

      if (!clusterValidPoissonNames(object@component@modelName))
      {stop("Invalid Poisson model name.")}
      return(TRUE)
    }
)

#' Initialize an instance of a MixAll S4 class.
#'
#' Initialization method of the [\code{\linkS4class{ClusterPoisson}}] class.
#' Used internally in the 'MixAll' package.
#'
#' @rdname initialize-methods
#' @keywords internal
setMethod(
    f="initialize",
    signature=c("ClusterPoisson"),
    definition=function(.Object, data=matrix(nrow=0, ncol=0), nbCluster=2, modelName="poisson_pk_ljk")
    {
      # for data
      if(missing(data)) {stop("data is mandatory in ClusterPoisson.")}
      .Object@component = new("ClusterPoissonComponent", data, nbCluster, modelName);
      .Object <- callNextMethod(.Object, nrow(.Object@component@data), nbCluster);
      # validate
      validObject(.Object)
      return(.Object)
    }
)

#' @rdname print-methods
#' @aliases print print,ClusterPoisson-method
#'
setMethod(
  f="print",
  signature=c("ClusterPoisson"),
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
#' @aliases show-ClusterPoisson,ClusterPoisson,ClusterPoisson-method
setMethod(
  f="show",
  signature=c("ClusterPoisson"),
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
  signature=c("ClusterPoisson"),
  function(object, ...)
  {
    cat("**************************************************************\n")
    callNextMethod()
    summary(object@component);
    cat("**************************************************************\n")
  }
)

#' Plotting of a class [\code{\linkS4class{ClusterPoisson}}]
#'
#' Plotting data from a [\code{\linkS4class{ClusterPoisson}}] object
#' using the estimated parameters and partition.
#'
#' @param x an object of class [\code{\linkS4class{ClusterPoisson}}]
#' @param y a list of variables to plot (subset). Variables names or indices.
#' If missingValues all the variables are represented.
#' @param ... further arguments passed to or from other methods
#'
#' @importFrom graphics plot
#' @aliases plot-ClusterPoisson, ClusterPoisson-method
#' @docType methods
#' @rdname plot-ClusterPoisson-method
#' @export
#'
#' @seealso \code{\link{plot}}
#' @examples
#'   ## for quantitative case
#' \dontrun{
#'   data(DebTrivedi)
#'   dt <- DebTrivedi[, c(1, 6,8, 15)]
#'   model <- clusterPoisson(iris[1:4], 3, strategy = clusterFastStrategy())
#'   plot(model)
#'   plot(model, c(1,2))
#'   }
#'
setMethod(
    f="plot",
    signature=c("ClusterPoisson"),
    function(x, y, ...)
    { # use generic plot
      .clusterPoissonPlot(x, y, .dPoisson,...);
    }
)

# wrapper of dnorm
# x a vector with the point
.dPoisson <- function(x, j, k, model)
{ dpois(x, lambda = (model@component@lambda)[k, j]) }

.clusterPoissonPlot <- function(model, y, ddensity,...)
{
  # total number of variable in the data set
  nbVariable = ncol(model@component@data);
  # no y => display all variables
  if (missing(y)) { y=1:nbVariable; }
  else # perform some check
  {
    if (is.numeric(y)) # numbers of the columns to plot are given
    {
      if (max(y)>nbVariable)
        stop("In .clusterPoissonPlot, y indices mismatch the data dimension")
    }
    else # names of the variables to plot are given
    {
      if ( sum(y %in% colnames(model@component@data))!= length(y) )
      { stop(cat("In .clusterPoissonPlot, unknown variable: ", paste(y[which(!(y %in% colnames(model@component@data)))]),"\n"))}
    }
  }
  # get old par
  op <- par(no.readonly = TRUE) # the whole list of settable par's.
  # cluster parameters
  par(mar = rep(2.5,4), cex = .75, oma = c(0, 0, 3, 0))        # margin and font size
  nbCol = length(y)                       # size of the matrix screen
  split.screen(c(nbCol, nbCol))           # create layout matrix screens
  col = model@zi+2;                       # color for each group
  pch = rep(1, length.out = length(col)); # circles
  pch[model@component@missing[,1]] = 3;            # + for missing values
  # create histograms on the diagonal
  for ( i in 1:nbCol )
  {
    screen(i+((i-1)*nbCol))   # sreen(i,i)
    xValues<-(min(model@component@data[,y[i]]):max(model@component@data[,y[i]]))
    density<-matrix(nrow=model@nbCluster, ncol=length(xValues))
    # loop over the clusters to generate densities
    for( k in 1:model@nbCluster )
    {  density[k,]<- ddensity(xValues, y[i], k, model);}
    # generate mixture density
    mixture<-apply(density,2,sum)
    if (is.numeric(y)) { xlab=colnames(model@component@data)[y[i]];}
    else               { xlab= y[i];}
    # TODO: check if xlab is empty
    main=paste("Histogram of",xlab)
    h<-hist(model@component@data[,y[i]], xlab=xlab, main=main, ...)
    # add on the histogram the estimated densities
    ratio<-max(h$counts)/max(mixture)
    density<-density*ratio
    mixture<-mixture*ratio
    lines(xValues,mixture,col="azure4", lty=1, lwd=4)
    for( k in 1:model@nbCluster )
    { lines(xValues, density[k,], col=k+1, lty=2, lwd=2)}
  }
  # add biplots
  if (nbCol>1)
  {
    for ( i in 2:nbCol )
    {
      if (is.numeric(y)) { xlab=colnames(model@component@data)[y[i]];}
      else               { xlab= y[i];}
      for( j in 1:(i-1) )
      {
        screen(j+((i-1)*nbCol)) # screen(i,j)
        if (is.numeric(y)) {ylab=colnames(model@component@data)[y[j]];}
        else {ylab= y[j];}
        plot(model@component@data[,y[j]], model@component@data[,y[i]], col=col, pch=pch, xlab=xlab, ylab=ylab, ...)
      }
    }
  }
#  mtext("Visualisation using latent logistic representation", outer = TRUE, cex = 1.5)
  close.screen(all.screens = TRUE)
  # restore plotting parameters
  par(op)
  invisible()
}
