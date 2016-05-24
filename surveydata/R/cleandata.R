# Functions to perform data cleanup


#
#  surveydata/R/cleandata.R by Andrie de Vries  Copyright (C) 2011-2012
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


#' Tests whether levels contain "Don't know".
#' 
#' Returns TRUE if x contains any instances of dk  
#' 
#' @param x Character vector or Factor 
#' @param dk Character vector, containing search terms, e.g. c("Don't know", "Don't Know")
#' @return TRUE or FALSE
#' @export 
#' @family Functions to clean data
#' @keywords "clean data"
hasDK <- function(x, dk="Don't Know"){
	ifelse(is.factor(x), 
			l <- levels(x),
			l <- unique(x)
	)
	any(l %in% dk)
}

#' Removes "Don't know" from levels and replaces with NA.
#' 
#' Tests the levels of x contain any instances of "Don't know".  If so, replaces these levels with NA
#' 
#' @param x Vector or Factor 
#' @param dk Character vector, containing search terms, e.g. c("Don't know", "Don't Know")
#' @return A factor with "Dont know" removed
#' @export 
#' @family Functions to clean data
#' @keywords "clean data"
removeDK <- function(x, dk="Don't Know"){
	if (hasDK(x, dk)){
    if(is.factor(x)){
  		l <- levels(x)
  		l[which(levels(x) %in% dk)] <- NA
  		x <- factor(x, levels=l)
    } else {
      pattern <- paste("^(", paste(dk, collapse="|"), ").?$", sep="")
      x <- gsub(pattern, "", x)
    }
  }  
  x
}

#' Removes "Do not know" and other similar words from factor levels in data frame.
#' 
#' Removes "Do not know" and other similar words from factor levels in data frame
#' 
#' @param x List or data frame 
#' @param dk Character vector, containing search terms, e.g. c("Do not know", "DK").  These terms will be replaced by NA. If NULL, will default to c("I don't know", "Don't Know", "Don't know", "Dont know" , "DK")
#' @seealso \code{\link{hasDK}} and \code{\link{removeDK}}
#' @return A data frame
#' @export
#' @family Functions to clean data
#' @keywords "clean data"
removeAllDK <- function(x, dk=NULL){
	if (is.null(dk)) dk <- c("I don't know", "Don't Know", "Don't know", "Dont know" , "DK")		
	newx <- lapply(x, removeDK, dk)
	n1 <- sum(as.numeric(lapply(x, hasDK, dk)))
	n2 <- sum(as.numeric(lapply(newx, hasDK, dk)))
	dk <- paste(dk, collapse=", ")
	message(paste("Removed", n1-n2,"instances of levels that equal [", dk, "]"))
	ret <- quickdf(newx)
	attributes(ret) <- attributes(x)
  class(ret) <- class(x)
	ret
}	

#' Fix level formatting of all question with Yes/No type answers.
#' 
#' @param x Surveyor data object
#' @export
#' @family Functions to clean data
#' @keywords "clean data"
leveltestSPSS <- function(x){
  ret <- FALSE
  if(inherits(x, "numeric")){
    if(!is.null(attributes(x)$value.labels)){
      if(all(attributes(x)$value.labels==c(1, 0))){
        ret <- TRUE
      }}}
  ret
}

#' Fix level formatting of all question with Yes/No type answers.
#' 
#' @param dat Surveyor data object
#' @export
#' @family Functions to clean data
#' @keywords "clean data"
fixLevels01SPSS <- function(dat){
  ret <- lapply(dat, function(x){
        if(leveltestSPSS(x)){
          x <- factor(x)
          levels(x) <- c("No", "Yes")
          x
        } else {
          x
        }
      }
  )
  ret <- plyr::quickdf(ret)
  attributes(ret)$variable.labels <- varlabels(dat)
  ret
}

#' Fix level formatting of all question with Yes/No type answers.
#' 
#' @param x Surveyor data object
#' @export
#' @keywords "clean data"
#' @family Functions to clean data
leveltestR <- function(x){
  ret <- FALSE
  if(inherits(x, "factor")){
    if(length(levels(x))==2){
      if(all(levels(x)==c("Yes", "Not selected"))){
        ret <- TRUE
      }}}
  ret
}

#' Fix level formatting of all question with Yes/No type answers.
#' 
#' @param dat surveydata object
#' @export
#' @family Functions to clean data
#' @keywords "clean data"
fixLevels01R <- function(dat){
  stopifnot(is.surveydata(dat))
  ret <- lapply(dat, function(x){
        if(leveltestR(x)){
          levels(x) <- c("Yes", "No")
          x
        } else {
          x
        }
      }
  )
  ret <- plyr::quickdf(ret)
  pattern(ret) <- pattern(dat)
  varlabels(ret) <- varlabels(dat)
  as.surveydata(ret)
}

#' Fix level formatting of all question with Yes/No type answers.
#' 
#' @param dat surveydata object
#' @param origin Either "R" or "SPSS"
#' @export
#' @family Functions to clean data
#' @keywords "clean data"
fixLevels01 <- function(dat, origin=c("R", "SPSS")){
  origin <- match.arg(origin)
  switch(origin,
      "R" = fixLevels01R(dat),
      "SPSS" = fixLevels01SPSS(dat))
}

#' Changes vector to ordered factor, adding NA levels if applicable.
#' 
#' @param x Vector
#' @export
#' @family Tools
qOrder <- function(x){
  #factor(x, level=levels(x), labels=levels(x), ordered=TRUE)
  if(any(is.na(x))){
    addNA(ordered(x))
  } else {
    ordered(x)
  }
}

#' Applies function only to named elements of a list.
#' 
#' This is useful to clean only some columns in a list (or data.frame or surveydata object). This is just a wrapper around \code{\link{lapply}} where only the named elements are changed.
#' @param x List
#' @param names Character vector identifying which elements of the list to apply FUN
#' @param FUN The function to apply.
#' @param ... Additional arguments passed to FUN
#' @export 
#' @family Tools
lapplyNames <- function(x, names, FUN, ...){
  oldClass <- class(x)
  index <- match(names, names(x))
  if(any(is.na(index))) stop(paste("Names not found:", paste(names[is.na(index)], collapse=", ")))
  x <- unclass(x)
  x[index] <- lapply(x[index], FUN, ...)
  class(x) <- oldClass
  x
}

#lapplyNames <- function(x, names, FUN, ...){
#  index <- match(names, names(x))
#  if(any(is.na(index))) stop(paste("Names not found:", paste(names[is.na(index)], collapse=", ")))
#  x[index, drop=TRUE] <- lapply(x[index, drop=TRUE], FUN, ...)
#  x
#}
