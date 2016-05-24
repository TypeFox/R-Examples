#' @name pts
#' @title The total phosphorus concentration and river discharge data of West Fork Cedar River at Finchford, Iowa.
#'
#' @description This data set contains the monthly total phosphorus concentration (P) and river discharge
#' (Q) of
#' West Fork Cedar River at Finchford, Iowa, USA, from 10/1998 to 10/2013. The P
#' data were collected under the ambient water quality program conducted by the
#' Iowa Department of Natural Resources (Iowa DNR), courtesy of Dr. K. E. Schilling from Iowa Geological Survey, 
#' University of Iowa. The Q data were obtained from the website of U.S.
#' Geological Survey. A gap from 09/2008 to 03/2009 in the data are due to
#' program suspension owing to lack 
#' of funding. The data contains \code{logP} and \code{logQ} which are the
#' logarithm of the original P and Q respectively. The P data are left censored,
#' with -1 or 0 in \code{ci} indicating that the corresponding observation is left censored or
#' observed. The lower censoring limits are stored in \code{lcl}. \code{season} is the quarter to which a month belong. 
#' @author Chao Wang (chao-wang@@uiowa.edu), 08/2015
NULL
