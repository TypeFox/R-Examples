## Experimental research in evolutionary computation
## author: thomas.bartz-beielstein@fh-koeln.de
## http://www.springer.com/3-540-32026-1
##
## Copyright (C) 2006-2010 T. Bartz-Beielstein and C. Lasarczyk
## This program is free software;
## you can redistribute it and/or modify it under the terms of the 
## GNU General Public License as published by the Free Software Foundation; 
## either version 3 of the License,
## or (at your option) any later version.
## This program is distributed in the hope that it will be useful, 
## but WITHOUT ANY WARRANTY; without even the implied warranty of 
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## See the GNU General Public License for more details.
## You should have received a copy of the GNU General Public License along
##  with this program; if not, see <http://www.gnu.org/licenses/>.
##

###################################################################################################
#' Get Repeats Function
#'
#' Helper Function for Spot's stats, evaluating the number of repeats
#' for a given parameter configuration
#'
#' @param configNumber one configuration number
#' @param configList list of configuration numbers 
#' @return number \code{repList} \cr
#' -  \code{repList} contains the number of repeats of the given config
#' @export
#' @keywords internal
###################################################################################################

spotGetRepeats <- function(configNumber, configList){
	listElements <- unique(configList)
	splitList <- (split(configList,factor(configList)))
	repList <- data.frame(cbind(listElements,reps=sapply(splitList, length)))
	repList[configNumber==listElements,2]
}


###################################################################################################
#' Get ALL Repeats Function
#'
#' Gets repeats of all configs in the list
#' 
#' @param configList config list
#' @return list \code{all} \cr
#' - \code{all} contains all repeats of the given config list
#' @export
#' @keywords internal
###################################################################################################
spotGetAllRepeats <- function(configList){
	all<-NULL;
	for (i in 1:length(configList)){
		all<-c(all,spotGetRepeats(configList[i],configList))
	}
	all
}

###################################################################################################
#' Helper function for transformations/normalizing
#'
#' Transforms natural to coded variables, used by spotPredictLm
#'
#' @param x parameter vector
#' @return vector \code{y} \cr
#' - \code{y} is the transformed/normalized version of \code{x}
#' @export
#' @keywords internal
###################################################################################################
spotHlpF.norm <- function(x){
	2* ( x - mean(x) ) / ( max(x) - min(x) )
}

###################################################################################################
#' Helper function for transformations/normalizing
#'
#' Transforms natural to coded variables
#'
#' @param x parameter vector
#' @param z vector (zmin, zmax)
#' @return vector \code{y} \cr
#' - \code{y} is the transformed/normalized version of \code{x}
#' @export
#' @keywords internal
###################################################################################################
spotHlpF.norminv <- function(x,z){
	x/2 * ( max(z) - min(z) ) +  mean(z) 
}