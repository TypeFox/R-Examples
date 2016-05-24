# Defines merge method

#
#  surveydata/R/merge.R by Andrie de Vries  Copyright (C) 2011-2012
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


#' Merges variable.labels attribute from two surveydata objects
#' 
#' Merges variable labels from two data objects.  The labels from dat1 takes precedence.
#' 
#' @param dat1 surveydata object
#' @param dat2 surveydata object
#' @param new_names A vector with names of the merged varlabels.  Defaults to the union of names of dat1 and dat2
#' @keywords internal
merge_varlabels <- function(dat1, dat2, new_names=union(names(dat1), names(dat2))){
  labels1 <- varlabels(dat1)
  labels2 <- varlabels(dat2)
  names(labels1) <- names(dat1)
  names(labels2) <- names(dat2)
  #merge(labels1, labels2)
  ret <- new_names
  names(ret) <- ret
  ret[names(labels2)] <- labels2
  ret[names(labels1)] <- labels1
  ret
}


#' Merge surveydata objects.
#' 
#' The base R merge will merge data but not all of the attributes.  This function also merges the variable.labels attribute.
#'
#' @name merge
#' @aliases merge merge.surveydata
#' @param x surveydata object
#' @param y surveydata object
#' @param ... Other parameters passed to \code{\link{merge}}
#' @method merge surveydata
#' @export
merge.surveydata <- function(x, y, ...){
  tmp <- merge(as.data.frame(x), as.data.frame(y), ...)
  newlabels <- merge_varlabels(x, y, new_names=names(tmp))
  varlabels(tmp) <- newlabels
  if(!identical(pattern(x), pattern(y))) warning("In merge of surveydata objects, patterns of objects differ")
  as.surveydata(tmp, ptn=pattern(x))
}


#' Combines surveydata object by columns.
#' 
#' @param ... surveydata objects
#' @param deparse.level ignored
#' @method cbind surveydata
#' @export
cbind.surveydata <- function(..., deparse.level=1){
  ptn <- pattern(..1)
  varlab <- do.call(c, lapply(list(...), varlabels))
  ret <- do.call(cbind.data.frame, list(...))
  varlabels(ret) <- varlab
  as.surveydata(ret)
}




