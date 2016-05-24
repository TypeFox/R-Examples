###  String functions

#
#  surveydata/R/strings.R by Andrie de Vries  Copyright (C) 2011-2012
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


#' Finds the common and unique elements in a character vector.
#' 
#' Function takes a character string as input and find the common and
#' unique elements.  Assumes that the common element is at start of string
#' 
#' @param string Character vector
#' @return list of common and unique strings 
#' @keywords string
#' @export 
#' @family Strings
#' @examples
#' test <- c("Q_1", "Q_2", "Q_3") 
#' strCommonUnique(test)$common
#' strCommonUnique(test)$unique
strCommonUnique <- function(string){
	x <- as.character(string)
	y <- string
	
	## Handles case with a single string element
	if (length(x) <= 1){
		return(list(common=x[1], unique=""))
	} 

	## Handles case where all elements are identical
	all_identical <- all(as.logical(lapply(x, function(f)x[1]==f)))
	if (all_identical){
		return(list(common=x[1], unique=rep("", length(x))))			
	}

	## Handles case where shortest element has length 0
	if (min(nchar(x))==0){
		return(list(common="", unique=x))
	}
	
	## Handles case where shortest element has length 1
	if (min(nchar(x))==1){
		x1 <- sapply(x, function(f){unlist(strsplit(f, NULL))[1]})
		all_identical <- all(as.logical(lapply(x1, function(f)x1[1]==f)))
		if (all_identical){
			return(
					list(common=substr(x[1], 1, 1), unique=substr(x, 2, nchar(x)))
			)			
		} else {
			return(
					list(common="", unique=x)
			)
		}	
	}
	
	
	# Make all strings the same length as shortest string
	x1 <- substr(x, 1, min(nchar(x)))
	# Create matrix of characters
	split <- lapply(x1, function(f){unlist(strsplit(f, NULL))})
	# Test which characters are identical
	identical <- sapply(split, function(f){f==split[[1]]}) ### aaply
	common <- apply(identical, 2, function(f){which(f==FALSE)[1]})  ### aaply
	mincommon <- min(common, na.rm=TRUE)-1
  #browser()
  if (mincommon <1){
		return(list(common="", unique=x))
	} else {
		return(list(
						common=substr(x[1], 1, mincommon),
						unique=substr(x, mincommon+1, nchar(x))
				))
	}
}




