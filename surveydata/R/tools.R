#
#  surveydata/R/tools.R by Andrie de Vries  Copyright (C) 2011-2012
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


#' Calculates at which questions repondents drop out.
#' 
#' The number of respondents for each question is calculated as the length of the vector, after omitting NA values.
#' 
#' @param x surveydata object, list or data.frame
#' @param summary If TRUE, returns a shortened vector that contains only the points where respondents drop out. Otherwise, returns the number of respondents for each question.
#' @return Named numeric vector of respondent counts
#' @export 
#' @examples
#' dropout(membersurvey[-(127:128)])
dropout <- function(x, summary=TRUE){
  len <- sapply(x, function(xx)length(na.omit(xx)))
  ll <- rev(cummax(rev(len)))
  len[c(1, 1+which(diff(ll) != 0))]
}

