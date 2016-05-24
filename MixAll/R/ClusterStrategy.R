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
#' @include ClusterInit.R
NULL

#-----------------------------------------------------------------------
#' Create an instance of [\code{\linkS4class{ClusterStrategy}}] class
#'
#' A strategy is a multistage empirical process for finding a
#' good estimate in the clustering estimation process.
#'
#' A strategy is a way to find a good estimate of the parameters of a mixture model
#' when using an EM algorithm or its variants. A ``try'' is composed of three stages
#' \itemize{
#'   \item \code{nbShortRun} short iterations of the initialization step and
#'    of the \code{EM}, \code{CEM}, \code{SEM} or \code{SemiSEM} algorithm.
#'   \item \code{nbInit} initializations using the [\code{\link{clusterInit}}]
#'   method.
#'   \item A long run of the \code{EM}, \code{CEM}, \code{SEM} or \code{SemiSEM} algorithm.
#' }
#' For example if \code{nbInit} is 5 and \code{nbShortRun} is also 5, there will
#' be 5 packets of 5 models initialized. In each packet, the best model will be
#' ameliorated using a short run. Among the 5 models ameliorated one will be
#' estimated until convergence using a long run. In total there were 25 initializations.
#'
#' The whole process can be repeated at least \code{nbTry} times. If a try
#' success, the estimated model is returned, otherwise an empty model is returned.
#'
#' @param nbInit Integer defining the number of initialization to try. Default value: 3.
#' @param initMethod Character string with the initialization method, see [\code{\link{clusterInit}}]$
#' for possible values. Default is "class".
#' @param initAlgo Character string with the algorithm to use in the initialization stage,
#' [\code{\link{clusterAlgo}}] for possible values. Default value: "EM".
#' @param nbInitIteration Integer defining the maximal number of iterations in initialization algorithm
#' if \code{initAlgo} = "EM", "CEM" or "SemiSEM". This is the number of iterations if \code{initAlgo} = "SEM".
#' Default value: 20.
#' @param initEpsilon Real defining the epsilon value for the algorithm.
#' \code{initEpsilon} is not used by the \code{SEM} algorithm. Default value: 0.01.
#'
#' @param nbShortRun Integer defining the number of short run to try
#' (the strategy launch an initialization before each short run). Default value: 5.
#' @param shortRunAlgo A character string with the algorithm to use in the short run stage
#' Default value: "EM".
#' @param nbShortIteration Integer defining the maximal number of iterations in the short runs
#' if \code{shortRunAlgo} = "EM", "CEM" or "semiSEM", or the number of iterations if \code{shortRunAlgo} = "SEM".
#' Default value: 100.
#' @param shortEpsilon Real defining the epsilon value for the algorithm.
#' epsilon is not used by the \code{SEM} algorithm. Default value: 1e-04.
#'
#' @param longRunAlgo A character string with the algorithm to use in the long run stage
#' Default value: "EM".
#' @param nbLongIteration  Integer defining the maximal number of iterations in the short runs
#' if \code{shortRunAlgo} = "EM", "CEM" or "SemiSEM", or the number of iterations if \code{shortRunAlgo} = "SEM".
#' Default value: 1000.
#' @param longEpsilon Real defining the epsilon value for the algorithm.
#' epsilon is not used by the \code{SEM} algorithm. Default value: 1e-07.
#'
#' @param nbTry number of estimation to attempt.
#'
#' @examples
#'    clusterStrategy()
#'    clusterStrategy(longRunAlgo= "CEM", nbLongIteration=100)
#'    clusterStrategy(nbTry = 1, nbInit= 1, shortRunAlgo= "SEM", nbShortIteration=100)
#'
#' @return a [\code{\linkS4class{ClusterStrategy}}] object
#' @author Serge Iovleff
#' @export
clusterStrategy <- function( nbTry =1
                           , nbInit= 5, initMethod="class", initAlgo= "EM", nbInitIteration=20, initEpsilon=0.01
                           , nbShortRun= 5, shortRunAlgo= "EM", nbShortIteration=100, shortEpsilon=1e-04
                           , longRunAlgo= "EM", nbLongIteration=1000, longEpsilon=1e-07
                           )
{
  # create init
  init = clusterInit(initMethod, nbInit, initAlgo, nbInitIteration, initEpsilon);
  # create shortAlgo
  shortAlgo = clusterAlgo(shortRunAlgo, nbShortIteration, shortEpsilon);
  # create longAlgo
  longAlgo = clusterAlgo(longRunAlgo, nbLongIteration, longEpsilon);
  # create strategy
  new("ClusterStrategy", nbTry =nbTry, nbShortRun =nbShortRun, initMethod =init, shortAlgo =shortAlgo, longAlgo =longAlgo);
}

#' @description
#' \code{clusterSemiSEMStrategy()} create an instance of [\code{\linkS4class{ClusterStrategy}}]
#' for users with many missing values.
#' @examples
#'    clusterSemiSEMStrategy()
#'
#' @rdname clusterStrategy
#' @export
clusterSemiSEMStrategy <- function()
{
  # create init
  initMethod = clusterInit("class", 5, "SemiSEM", 20, 0);
  # create shortAlgo
  shortAlgo = clusterAlgo("SemiSEM", 50, 0);
  # create longAlgo
  longAlgo = clusterAlgo("SemiSEM", 400, 0);
  # create strategy
  new("ClusterStrategy", nbTry= 2, nbShortRun= 5, initMethod= initMethod, shortAlgo= shortAlgo, longAlgo= longAlgo);
}

#' @description
#' \code{clusterFastStrategy()} create an instance of [\code{\linkS4class{ClusterStrategy}}] for impatient user.
#' @examples
#'    clusterFastStrategy()
#'
#' @rdname clusterStrategy
#' @export
clusterFastStrategy <- function()
{
  # create init
  initMethod = clusterInit("class", 3, "EM", 5, 0.01);
  # create shortAlgo
  shortAlgo = clusterAlgo("CEM", 10, 0.001);
  # create longAlgo
  longAlgo = clusterAlgo("EM", 100, 1e-07);
  # create strategy
  new("ClusterStrategy", nbTry= 1, nbShortRun= 2, initMethod= initMethod, shortAlgo= shortAlgo, longAlgo= longAlgo);
}

#-----------------------------------------------------------------------
#' Constructor of [\code{\linkS4class{ClusterStrategy}}] class
#'
#' This class encapsulate the parameters of the clustering estimation strategies.
#'
#'   @slot nbTry Integer defining the number of tries. Default value: 1.
#'   @slot nbShortRun Integer defining the number of short run. Recall that the
#'   strategy launch an initialization before each short run. Default value is 5.
#'   @slot initMethod A [\code{\linkS4class{ClusterInit}}] object defining the way to
#'   initialize the estimation method. Default value is [\code{\link{ClusterInit}}].
#'   @slot shortAlgo A [\code{\linkS4class{ClusterAlgo}}] object defining the algorithm
#'   to use during the short runs of the estimation method. Default value is
#'   \code{clusterAlgo("EM",100,1e-04)}.
#'   @slot longAlgo A [\code{\linkS4class{ClusterAlgo}}] object defining the algorithm
#'   to use during the long run of the estimation method.  Default value is
#'   \code{clusterAlgo("EM",1000,1e-07)}.
#'
#' @examples
#'   new("ClusterStrategy")
#'   shortAlgo=clusterAlgo("SEM",1000)
#'   longAlgo =clusterAlgo("SemiSEM",200,1e-07)
#'   new("ClusterStrategy", shortAlgo=shortAlgo, longAlgo=longAlgo)
#'   getSlots("ClusterStrategy")
#'
#' @author Serge Iovleff
#'
#' @name ClusterStrategy
#' @rdname ClusterStrategy-class
#' @aliases ClusterStrategy-class
#' @exportClass ClusterStrategy
setClass(
    Class="ClusterStrategy",
    slots=c( nbTry = "numeric", nbShortRun = "numeric", initMethod = "ClusterInit", shortAlgo="ClusterAlgo", longAlgo="ClusterAlgo" ),
    prototype=list(  nbTry = 1, nbShortRun=5, initMethod = clusterInit(), shortAlgo=clusterAlgo("EM",100,1e-04), longAlgo=clusterAlgo("EM",1000,1e-07)),
    # validity function
    validity=function(object)
    {
      if (round(object@nbTry)!=object@nbTry)
      {stop("nbTry must be an integer.")}
      if( object@nbTry < 1 ) # can't be zero
      {stop("nbTry must be strictly greater than 0");}
      if (round(object@nbShortRun)!=object@nbShortRun)
      {stop("nbShortRun must be an integer.")}
      if(class(object@initMethod)[1] != "ClusterInit")
      {stop("initMethod is not of a Cluster Initialization method (must be an instance of the class ClusterInit).")}
      if(class(object@shortAlgo)[1] != "ClusterAlgo")
      {stop("shortAlgo is not of a Cluster algorithm (must be an instance of the class ClusterAlgo).")}
      if(class(object@longAlgo)[1] != "ClusterAlgo")
      {stop("longAlgo is not of a Cluster algorithm (must be an instance of the class ClusterAlgo).")}
      return(TRUE)
    }
)

#' Initialize an instance of a MixAll S4 class.
#'
#' Initialization method of the [\code{\linkS4class{ClusterStrategy}}] class.
#' Used internally in the 'MixAll' package.
#'
#' @keywords internal
#' @rdname initialize-methods
#'
setMethod(
  f="initialize",
  signature=c("ClusterStrategy"),
  definition=function(.Object, nbTry, nbShortRun, initMethod, shortAlgo, longAlgo)
  {
    # for nbtry
    if(missing(nbTry)) {.Object@nbTry<-1}
    else  {.Object@nbTry<-nbTry}
    # for nbtry
    if(missing(nbShortRun)) {.Object@nbShortRun<-5}
    else  {.Object@nbShortRun<-nbShortRun}
    # for initMethod
    if( missing(initMethod) ){ .Object@initMethod<-clusterInit() }
    else{.Object@initMethod<-initMethod}
    # for shortAlgo
    if(missing(shortAlgo)){ .Object@shortAlgo<-clusterAlgo("EM", 100, 1e-04) }
    else{.Object@shortAlgo<-shortAlgo}
    # for longAlgo
    if(missing(longAlgo)){ .Object@longAlgo<-clusterAlgo("EM", 1000, 1e-07) }
    else{.Object@longAlgo<-longAlgo}
    # validate
    validObject(.Object)
    return(.Object)
  }
)

#' @aliases print print-strategy,ClusterStrategy-method
#' @rdname print-methods
setMethod(
  f="print",
  signature=c("ClusterStrategy"),
  function(x,...){
    cat("****************************************\n")
    cat("*** Cluster Strategy:\n")
    cat("* number of try         = ", x@nbTry, "\n")
    cat("* number of short run   = ", x@nbShortRun, "\n")
    cat("****************************************\n")
    cat("*** Initialization :\n")
    cat("* method = ", x@initMethod@method, "\n")
    cat("* number of init       = ", x@initMethod@nbInit, "\n")
    cat("* algorithm            = ", x@initMethod@algo@algo, "\n")
    cat("* number of iterations = ", x@initMethod@algo@nbIteration, "\n")
    cat("* epsilon              = ", x@initMethod@algo@epsilon, "\n")
    cat("****************************************\n")
    cat("*** short algorithm :\n")
    cat("* algorithm            = ", x@shortAlgo@algo, "\n")
    cat("* number of iterations = ", x@shortAlgo@nbIteration, "\n")
    cat("* epsilon              = ", x@shortAlgo@epsilon, "\n")
    cat("****************************************\n")
    cat("*** long algorithm :\n")
    cat("* algorithm            = ", x@longAlgo@algo, "\n")
    cat("* number of iterations = ", x@longAlgo@nbIteration, "\n")
    cat("* epsilon              = ", x@longAlgo@epsilon, "\n")
    cat("****************************************\n")
  }
)

#' @rdname show-methods
#' @aliases show show-strategy,ClusterStrategy-method
setMethod(
  f="show",
  signature=c("ClusterStrategy"),
  function(object){
    cat("****************************************\n")
    cat("*** Cluster Strategy:\n")
    cat("* number of try         = ", object@nbTry, "\n")
    cat("* number of short run   = ", object@nbShortRun, "\n")
    cat("****************************************\n")
    cat("*** Initialization :\n")
    cat("* method = ", object@initMethod@method, "\n")
    cat("* number of init       = ", object@initMethod@nbInit, "\n")
    cat("* algorithm            = ", object@initMethod@algo@algo, "\n")
    cat("* number of iterations = ", object@initMethod@algo@nbIteration, "\n")
    cat("* epsilon              = ", object@initMethod@algo@epsilon, "\n")
    cat("****************************************\n")
    cat("*** short Algorithm :\n")
    cat("* algorithm            = ", object@shortAlgo@algo, "\n")
    cat("* number of iterations = ", object@shortAlgo@nbIteration, "\n")
    cat("* epsilon              = ", object@shortAlgo@epsilon, "\n")
    cat("****************************************\n")
    cat("*** long algorithm :\n")
    cat("* algorithm            = ", object@longAlgo@algo, "\n")
    cat("* number of iterations = ", object@longAlgo@nbIteration, "\n")
    cat("* epsilon              = ", object@longAlgo@epsilon, "\n")
    cat("****************************************\n")
  }
)

#' @rdname extract-methods
#' @aliases [,ClusterStrategy-method
#'
setMethod(
  f="[",
  signature(x = "ClusterStrategy"),
  definition=function(x,i,j,drop){
    if ( missing(j) ){
      switch(EXPR=i,
        "nbTry"={return(x@nbTry)},
        "nbShortRun"={return(x@nbShortRun)},
        "initMethod"={return(x@initMethod)},
        "shortAlgo"={return(x@shortAlgo)},
        "longAlgo"={return(x@longAlgo)},
        stop("This attribute doesn't exist !")
      )
    }else{
      stop("This attribute is not a list !")
    }
  }
)

#-----------------------------------------------------------------------
#' @name [
# @docType methods
#' @rdname extract-methods
#' @aliases [<-,ClusterStrategy-method
setReplaceMethod(
  f="[",
  signature(x = "ClusterStrategy"),
  definition=function(x,i,j,value){
    if ( missing(j) )
    {
      switch(EXPR=i,
             "nbTry"={x@nbTry<-value},
             "nbShortRun"={x@nbShortRun<-value},
             "initMethod"={x@initMethod<-value},
             "shortAlgo"={x@shortAlgo<-value},
             "longAlgo"={x@longAlgo<-value},
             stop("This attribute doesn't exist !")
      )
    }else{
      stop("This attribute is not a list !")
    }
    validObject(x)
    return(x)
  }
)
