

## ## .. content for \description{} (no empty lines) ..
## ##
## ## .. content for \details{} ..
## ## @title scale.and.crossprod
## ## @param x 
## ## @examples
## ## x <- matrix(sample(1:12),3)
## ## scale.and.crossprod(x)
## ## cor(x)
## scale.and.crossprod <-
##   function(x)
## {
##   x <- .scale.data(x,scale="column")
##   ret <- crossprod(x)/2
##   ret[is.na(ret)] <- 0
##   ret
## }



is.dist <-
  function(d)
{
  inherits(d,"dist")
}



##' \code{gdist} computes and returns the distance matrix computed by using user-defined distance measure.
##'
##' \code{is.dist} tests if its argument is a `dist' object.
##'
##' The distance (or dissimilarity) function (\code{FUN}) can be any distance measure applied to \code{x}.
##' For instance, \code{"euclidean"}, \code{"maximum"}, \code{"manhattan"},\code{"canberra"},
##' \code{"binary"}, \code{"minkowski"}, "correlation.of.variables", "correlation.of.observations" or
##' \code{gmdm}. "correlation.of.variables" computes the correlation distance of
##' the variables (the columns); all the other compute the distances between
##' the observations (the rows) of a data matrix.
##' @title Generalized Distance Matrix Computation
##' @name gdist
##' @aliases gdist is.dist
##' @usage
##' gdist(x,method="euclidean",MoreArgs=NULL,diag=FALSE,upper=FALSE)
##'
##' is.dist(d)
##' 
##' @param x a numeric matrix, data frame or `dist' object.
##' @param method the distance measure to be used. This can either be one of
##' the methods used in \code{dist} (see \code{help("dist", package="stats")})
##' or \code{"correlation"}, \code{"correlation.of.observations"} and
##' \code{"correlation.of.variables"}. In addition, user-defined distance measure
##' are also allowed, which returns a \emph{dist} object and should at least
##' have attributes \emph{"Size"} and \emph{"Labels"}.
##' @param MoreArgs a list of other arguments to be passed to \code{gdist}.
##' @param diag logical value indicating whether the diagonal of the distance matrix should be 
##' printed by \code{print.dist}.
##' @param upper logical value indicating whether the upper triangle of the distance matrix should be
##' printed by \code{print.dist}.
##' @param d an R object.
##' @return
##' \code{gdist} returns an object of `dist'.\cr
##' \code{is.dist} returns a logical value whether an object is `dist'.\cr
##' @examples
##' ## load library
##' require("GMD")
##' require(cluster)
##' 
##' ## compute distance using Euclidean metric (default)
##' data(ruspini)
##' x <- gdist(ruspini)
##' 
##' ## see a dendrogram result by hierarchical clustering
##' dev.new(width=12, height=6)
##' plot(hclust(x),
##'      main="Cluster Dendrogram of Ruspini data",
##'      xlab="Observations")
##' 
##' ## convert to a distance matrix
##' m <- as.matrix(x)
##' 
##' ## convert from a distance matrix
##' d <- as.dist(m)
##' stopifnot(d == x)
##' 
##' ## Use correlations between variables "as distance"
##' data(USJudgeRatings)
##' dd <- gdist(x=USJudgeRatings,method="correlation.of.variables")
##' dev.new(width=12, height=6)
##' plot(hclust(dd),
##'      main="Cluster Dendrogram of USJudgeRatings data",
##'      xlab="Variables")
##' 
gdist <-
  function(x,
           method="euclidean",
           MoreArgs=NULL,
           diag=FALSE,
           upper=FALSE
           )
{
  if(method %in% c("correlation","correlation.of.observations")){
    FUN <- function(x,...){
      as.dist(1-cor(t(x),y=NULL,...),diag=diag,upper=upper)}
    if (.invalid(MoreArgs)) MoreArgs=list(method="pearson",use="everything")
  } else if(method %in% c("correlation.of.variables")){
    FUN <- function(x,...){
      as.dist(1-cor(x,y=NULL,...),diag=diag,upper=upper)}
    if (.invalid(MoreArgs)) MoreArgs=list(method="pearson",use="everything")
  }

  COMMON_METHODS <-
    c("euclidean","maximum",
      "manhattan","canberra",
      "binary","minkowski"
      )
  if(method %in% COMMON_METHODS){    
    d <- dist(x=x,method=method,diag=diag,upper=upper,p=MoreArgs$p)
  } else if (method %in% c("correlation","correlation.of.observations","correlation.of.variables")){
    ##d <- .call.FUN(FUN,x,MoreArgs)
    d <- FUN(x,method=MoreArgs$method,use=MoreArgs$use)
    attr(d,"method") <- method
  } else {
    FUN <- match.fun(method)
    MoreArgs[["diag"]] <- diag
    MoreArgs[["upper"]] <- upper
    d <- .call.FUN(FUN,x,MoreArgs)

    ## check attributes of the dist object ##
    if(is.null(attr(d,"method"))){
      attr(d,"method") <- method
    }
    if(is.null(attr(d,"call"))){
      attr(d,"call") <- match.call()
    }
    if(is.null(attr(d,"Size"))){
      warning(sprintf("The `dist' object returned by %s does not contain a specified attribute of `Size'.",attr(d,"method")))
    }
    if(is.null(attr(d,"Labels"))){
      warning(sprintf("The `dist' object returned by %s does not contain a specified attribute of `Labels'.",attr(d,"method")))
    }
  }
  attr(d,"Diag") <- diag
  attr(d,"Upper") <- upper
  class(d) <- "dist"
  return(d)
}


