# Subsetting code for surveydata objects

#
#  surveydata/R/extract.R by Andrie de Vries  Copyright (C) 2011-2012
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 or 3 of the License
#  (at your option).
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/
#


#' Extract or replace subsets of surveydata, ensuring that the varlabels stay in synch.
#'
#' The surveydata package makes it easy to extract specific questions from a surveydata object. Because survey data typically has question names like Q1_a, Q1_b, Q1_c the extract method for a surveydata object makes it easy to extract all columns by simply specifing "Q1" as the argument to the column index.
#' 
#' Extraction is similar to data frames, with three important exceptions:
#' \itemize{
#' \item{The column argument j is evaluated using \code{\link{which.q}} and will return all columns where the column names match the \code{\link{pattern}}.}
#' \item{The drop argument is FALSE. Thus the result will always be a surveydata object, even if only a single column is returned.}
#' \item{All extraction methods retain the \code{\link{pattern}} and \code{\link{varlabels}} arguments.}
#' }
#' 
#' @name Extract
#' @param i row index
#' @param j column index
#' @param drop logical. Passed to \code{\link{[.data.frame}}. Note that the default is FALSE.
#' @param ... Other arguments passed to \code{[.data.frame}
#' @export 
#' @aliases [ [.surveydata
#' @method [ surveydata
#' @example /inst/examples/example-extract.R
"[.surveydata" <- function(x, i, j, drop=FALSE){
  name <- NULL
  has.drop <- !missing(drop)
  Narg <- nargs() - (has.drop) -1
  has.i <- !missing(i)
  has.j <- !missing(j)
  
  if(!has.i & !has.j) return(x)
  
  if(Narg >= 1L & has.j){ 
    name <- j
    if (is.character(j)) {
      w <- which.q(x, j) 
      if(length(w)!=0) {
        j <- name <- w
      }
    } else {
      name <- j
    }
  } else { #!has.j
    name <- seq_along(x)
  }
  
  if(Narg==1L & has.i) {
    drop <- NULL
    name <- i
    if (is.character(i)) {
      w <- which.q(x, i)  
      if(length(w)!=0) {
        i <- name <- w
      } else {
        name <- i
      }
    } else { 
      name <- i
    } 
  } 
  
  if(is.null(drop)){
    ret <- NextMethod("[<-")
  } else {
    ret <- NextMethod("[<-", drop=drop)
  }
  
  if(is.data.frame(ret)){
    varlabels(ret) <- varlabels(x)[name]
    as.surveydata(ret, ptn=pattern(x), renameVarlabels=FALSE)
  } else {
    ret
  }
}



#' @rdname Extract 
#' @method [<- surveydata
#' @usage \method{[}{surveydata}(x, i, j) <- value
#' @export 
"[<-.surveydata" <- function(x, i, j, value){
  
  has.value <- !missing(value)
  Narg <- nargs() - (has.value) - 1
  
  has.i <- !missing(i)
  has.j <- !missing(j)
  if(Narg >= 1L & has.j){ 
    if(is.character(j)) {
      newname <- j
      w <- which.q(x, j)
      if(length(w)!=0) {
        j <- w
        name <- j
      } else {
        name <- newname
      }
      
    } else {
      name <- j
    }
  }
  
  if(Narg==1L & has.i) {
    if(is.character(i)) {
      newname <- i
      w <- which.q(x, i)
      if(length(w)!=0) {
        i <- w
        name <- i
      } else {
        name <- newname
      }
      
    } else {
      name <- i
    }
  }
  
  xorig <- x
  
  labels <- varlabels(x)
  if(is.null(value)){
    labels <- labels[-name]
  }
  if(length(w)==0){
    labels[newname] <- newname
  }  
  
  ret <- NextMethod("[<-")
  varlabels(ret) <- labels
  as.surveydata(ret, ptn=pattern(xorig), renameVarlabels=FALSE)
}


#' @rdname Extract 
#' @aliases $<- $<-.surveydata
#' @param x surveydata object
#' @param name Names of columns
#' @param value New value
#' @method $<- surveydata
#' @usage \method{$}{surveydata}(x, name) <- value
#' @export 
#' @seealso \code{\link{surveydata-package}}, \code{\link{varlabels}}
"$<-.surveydata" <- function(x, name, value){
  labels <- varlabels(x)
  if(is.null(value)){
    labels <- labels[names(labels)!=name]
  }
  if(length(grep(name, names(x)))==0){
    labels[name] <- name
  }  
  x <- as.data.frame(x)
  x <- NextMethod("$<-")
  varlabels(x) <- labels
  as.surveydata(x, renameVarlabels=FALSE)
}


