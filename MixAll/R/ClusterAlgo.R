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
#' Create an instance of the [\code{\linkS4class{ClusterAlgo}}] class
#'
#' There is three algorithms and two stopping rules possibles for an algorithm.
#' \itemize{
#'        \item Algorithms:
#'           \itemize{
#'               \item \code{EM}  {The Expectation Maximisation algorithm}
#'               \item \code{CEM} {The Classification EM algorithm}
#'               \item \code{SEM} {The Stochastic EM algorithm}
#'               \item \code{SemiSEM} {The Semi-Stochastic EM algorithm}
#'           }
#'        \item Stopping rules:
#'           \itemize{
#'               \item \code{nbIteration} {Set the maximum number of iterations.}
#'               \item \code{epsilon} {Set relative increase of the log-likelihood criterion.}
#'           }
#'        \item Default values are \eqn{200} \code{nbIteration} of \code{EM}
#'         with an \code{epsilon} value of \eqn{1.e-8}.
#' }
#' The epsilon value is not used when the algorithm is "SEM" or "SemiSEM".
#'
#' @param algo character string with the estimation algorithm.
#' Possible values are "EM", "SEM", "CEM", "SemiSEM". Default value is "EM".
#' @param nbIteration Integer defining the maximal number of iterations. Default value is 200.
#' @param epsilon Real defining the epsilon value for the algorithm. Not used
#'  by the "SEM" and "SemiSEM" algorithms. Default value is 1.e-7.
#'
#' @examples
#' clusterAlgo()
#' clusterAlgo(algo="SEM", nbIteration=50)
#' clusterAlgo(algo="CEM", epsilon = 1e-06)
#'
#' @return a [\code{\linkS4class{ClusterAlgo}}] object
#' @author Serge Iovleff
#' @export
#'
clusterAlgo <- function( algo="EM", nbIteration=200, epsilon=1e-07)
{
  # check algo
  if ( sum(algo %in% c("EM","SEM","CEM","SemiSEM")) != 1 )
  {  stop("algo is not valid. See ?clusterAlgo for the list of available algorithms.")}
  # check nbIteration
  if (!is.numeric(nbIteration))
  {stop("nbIteration must be an integer.")}
  if (round(nbIteration)!= nbIteration)
  {stop("nbIteration must be an integer.")}
  if( nbIteration < 0 ) # can be zero (no iterations)
  { stop("nbIteration must be positive or zero.")}
  # check epsilon
  if (algo != "SEM" && algo != "SemiSEM")
  {
    if (!is.double(epsilon) )
    {  stop("epsilon must be a scalar.")}
    if( epsilon < 0.)
    {  stop("epsilon must be positive.")}
  }
  new("ClusterAlgo", algo=algo, nbIteration=nbIteration, epsilon=epsilon)
}

#-----------------------------------------------------------------------
#' @title [\code{\linkS4class{ClusterAlgo}}] class for Cluster algorithms.
#'
#' @description
#' This class encapsulates the parameters of clustering estimation algorithms
#' methods.
#'
#' @slot algo A character string with the algorithm.
#' Possible values: "SEM", "CEM", "EM", "SemiSEM". Default value: "EM".
#' @slot nbIteration Integer defining the maximal number of iterations. Default value: 200.
#' @slot epsilon real defining the epsilon value for the algorithm. epsilon is
#' note used if \code{algo} is "SEM" or "SemiSEM". Default value: 1e-07.
#'
#' @examples
#' getSlots("ClusterAlgo")
#' new("ClusterAlgo")
#' new("ClusterAlgo", algo="SEM", nbIteration=1000)
#'
#' @name ClusterAlgo
#' @rdname ClusterAlgo-class
#' @aliases ClusterAlgo-class
#' @exportClass ClusterAlgo
setClass (
  Class= "ClusterAlgo",
  representation(algo="character", nbIteration = "numeric", epsilon = "numeric"),
  prototype=list(algo="EM", nbIteration = 200, epsilon = 1e-07),
  # validity function
  validity = function(object)
  {
    # for algo
    if ( sum(object@algo %in% c("EM","SEM","CEM","SemiSEM")) != 1 )
    {  stop("Algorithm is not valid. See ?ClusterAlgo for the list of available algorithms.")}
    # for nbIteration
    if (!is.numeric(object@nbIteration))
    {stop("nbIteration must be an integer.")}
    if (round(object@nbIteration)!=object@nbIteration)
    {stop("nbIteration must be an integer.")}
    if( object@nbIteration < 0 ) # can be zero (no iterations)
    {  stop("nbIteration must be positive or zero.")}
    # for epsilon
    if (object@algo != "SEM" && object@algo != "SemiSEM")
    {
      if (!is.double(object@epsilon) )
      {  stop("epsilon must be a scalar.")}
      if( object@epsilon < 0.)
      {  stop("epsilon must be positive.")}
    }
    return(TRUE)
  }
)

#' Initialize an instance of a MixAll S4 class.
#'
#' Initialization method of the [\code{\linkS4class{ClusterAlgo}}] class.
#' Used internally in the `MixAll' package.
#'
#' @keywords internal
#' @rdname initialize-methods
setMethod(
    f="initialize",
    signature=c("ClusterAlgo"),
    definition=function(.Object,algo,nbIteration,epsilon)
    {
      # for algo
      if(missing(algo)) {.Object@algo<-"EM"}
      else  {.Object@algo<-algo}
      # for epsilon
      if( missing(epsilon) ){ .Object@epsilon<-1e-07 }
      else{.Object@epsilon<-epsilon}
      # for nbIteration
      if(missing(nbIteration)){ .Object@nbIteration<-200 }
      else{.Object@nbIteration<-nbIteration}
      # validate
      validObject(.Object)
      return(.Object)
    }
)

#' @rdname print-methods
#' @aliases print print-algo,ClusterAlgo-method
setMethod(
    f="print",
    signature=c("ClusterAlgo"),
    function(x,...){
      cat("****************************************\n")
      cat("*** MixAll ClusterAlgo:\n")
      cat("* algorithm            = ", x@algo, "\n")
      cat("* number of iterations = ", x@nbIteration, "\n")
      cat("* epsilon              = ", x@epsilon, "\n")
      cat("****************************************\n")
    }
)

#' @rdname show-methods
#' @aliases show show-algo,ClusterAlgo-method
setMethod(
    f="show",
    signature=c("ClusterAlgo"),
    function(object){
      cat("****************************************\n")
      cat("*** MixAll ClusterAlgo:\n")
      cat("* algorithm            = ", object@algo, "\n")
      cat("* number of iterations = ", object@nbIteration, "\n")
      cat("* epsilon              = ", object@epsilon, "\n")
      cat("****************************************\n")
    }
)

#' @rdname extract-methods
#' @aliases [,ClusterAlgo-method
setMethod(
  f="[",
  signature(x = "ClusterAlgo"),
  definition=function(x, i, j, drop){
    if ( missing(j) ){
      switch(EXPR=i,
          "algo"={return(x@algo)},
          "nbIteration"={return(x@nbIteration)},
          "epsilon"={return(x@epsilon)},
          stop("This attribute doesn't exist !")
      )
    }else{
      stop("This attribute is not a list !")
    }
  }
)

#' @name [
# @docType methods
#' @rdname extract-methods
#' @aliases [<-,ClusterAlgo-method
setReplaceMethod(
    f="[",
    signature(x = "ClusterAlgo"),
    definition=function(x,i,j,value){
      if ( missing(j) )
      {
        switch(EXPR=i,
            "algo"={x@algo<-value},
            "nbIteration"={x@nbIteration<-value},
            "epsilon"={x@epsilon<-value},
            stop("This attribute doesn't exist !")
        )
      }else{
        stop("This attribute is not a list !")
      }
      validObject(x)
      return(x)
    }
)
