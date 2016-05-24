## define class for traits
ptypes <- c("multitype","binary","continuous","DNA","RNA","aacid",
            "other","unknown")

##' Class "pdata"
##'
##' Data class for phylo4d objects
##'
##'
##' @name pdata-class
##' @aliases ptypes pdata-class [<-,pdata-method [,pdata-method
##' [,pdata,ANY,ANY,ANY-method [[,pdata-method [[<-,pdata-method
##' [[,pdata,ANY,ANY-method [[,pdata,ANY,missing-method
##' @docType class
##' @section Objects from the Class: Objects can be created by calls of the form
##' \code{new("pdata", ...)}.
##' @author Ben Bolker
##' @keywords classes
####  @export
setClass("pdata", representation(data="data.frame",
                                 type="factor",
                                 comment="character",
                                 metadata="list"),
         prototype=list(data=data.frame(),type=factor(),
           comment=character(0),metadata=list()))

## pdata constructor


##' Constructor for pdata (phylogenetic data) class
##'
##' Combine data, type, comments, and metadata information to create a new pdata
##' object, or check such an object for consistency
##'
##'
##' @aliases pdata check_pdata
##' @param data a data frame
##' @param type a factor with levels as specified by \linkS4class{pdata}, the
##' same length as \code{ncol(data)}
##' @param comment a character vector, the same length as \code{ncol(data)}
##' @param metadata an arbitrary list
## @param object an object of class \code{pdata}
##' @return An object of class \code{pdata}
##' @author Ben Bolker
##' @seealso \linkS4class{pdata}
##' @keywords misc
pdata <- function(data,type,comment,metadata) {
  nvar <- ncol(data)
  if (missing(type)) {
    type <- factor(rep("unknown",nvar),levels=ptypes)
  }
  if (length(type)==1) type <- rep(type,length.out=nvar)
  type <- factor(as.character(type),levels=ptypes)
  if (length(comment)==1) comment <- rep(comment,length.out=nvar)
  obj <- new("pdata",data=data,type=type,comment=comment,metadata)
  check_pdata(obj)
  obj
}


check_pdata <- function(object) {
    nvar <- ncol(object@data)
    badlevels <- levels(object@type)[!levels(object@type) %in% ptypes]
    if (length(badlevels)>0)
      stop(paste("bad levels in types:",paste(badlevels,collapse=",")))
    if (length(object@comment)>1 && length(object@comment)!=nvar) {
        stop("wrong number of comments")
    }
    if (length(object@type)>1 && length(object@type)!=nvar) {
        stop("wrong number of types")
    }
}

## setMethod("[","pdata",function(x,i, j,...,drop=FALSE) {
##   xd <- x@data[i,j,...,drop=drop]
##   xd2 <- as.data.frame(xd)
##   xd2
## })

## #### @exportMethod [<-
## setGeneric("[<-")

## setMethod("[<-","pdata",function(x,i, j,...,drop=FALSE,value) {
##   "[<-"(x@data,i,j,...,drop=drop,value)
## })

## ### @exportMethod [[
## setGeneric("[[")
## setMethod("[[","pdata",
##           function(x,i,j,...,exact=NA) {
##             x@data[[i,j,...,exact=exact]]
##           })

## #### @exportMethod [[<-
## setGeneric("[[<-")
## setMethod("[[<-","pdata",
##           function(x,i,j,...,exact=NA,value) {
##             "[[<-"(x@data,i,j,...,exact=exact,value)
##           })

## setMethod("plot",signature(x="pdata",y="missing"), function(x,...){
##     return(plot(x@data, ...))
## }) # end plot phylo4


## od = data.frame(a=1:3,b=4:6)
## z = new("pdata",
##   data=od,type=factor("a","b"),
##   comment=c("",""),metadata=list())

## z[2,]
## z[,"a"]
## z[[2]]

## test conflict resolution error

#######
### old code retrieved from misc/ folder

## setClass("pdata", representation(x="vector", y="vector"))
## setMethod("[","pdata",function(x,i, j,...,drop=TRUE)new("pdata",x=x@x[i],y=x@y[i]))

# x <- new("pdata", x=c("a","b", "c", "d", "3"), y=c(1:5))
#>x[c(2,4)]
#An object of class pdata
#Slot "x":
#[1] "b" "d"
#
#Slot "y":
#[1] 2 4



# doesn't work
#setClass("track", representation("list", comment="character", metadata="vector"), contains="list", prototype(list(), comment="", metadata=NA))
#setMethod("[","track",function(x,i, j,...,drop=TRUE)new("track", list(lapply(x, function(x, i, j, ..., drop=TRUE) x@.Data[i]))))

# this works, how to incorporate into method above?
#> lapply(x, function(x, i=2, j, ..., drop=TRUE) x@.Data[i])
#$x
#[1] "b"

#$y
#[1] 2

# this works, but list structure is destroyed
#> mapply(function(x, i, j, ..., drop=TRUE) x@.Data[i], x, 2)
#  x   y
#"b" "2"
