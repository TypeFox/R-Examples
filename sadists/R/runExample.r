# Copyright 2014-2015 Steven E. Pav. All Rights Reserved.
# Author: Steven E. Pav
#
# This file is part of sadists.
#
# sadists is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# sadists is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with sadists.  If not, see <http://www.gnu.org/licenses/>.

# Created: 2015.03.16
# Copyright: Steven E. Pav, 2015
# Author: Steven E. Pav
# Comments: Steven E. Pav


#' @title Run Shiny Application
#'
#' @description 
#'
#' Runs a shiny application which draws from the given distributions, then
#' illustrates the fidelity of the density, CDF, and quantile functions.
#'
#' @details
#'
#' Launches shiny applications, and optionally, your system's web browser.
#' Draws are taken from the random variable, and d-d, q-q, and p-p plots
#' are available.
#'
#' @usage
#'
#' runExample(port=NULL,launch.browser=TRUE,
#'  host=getOption('shiny.host','127.0.0.1'),display.mode='auto')
#'
#' @inheritParams shiny::runApp
#' @references
#'
#' Attali, D. "Supplementing your R package with a shiny app."
#' \url{http://deanattali.com/2015/04/21/r-package-shiny-app/}
#'
#' @export 
#' @examples 
#' \dontrun{
#' runExample(launch.browser=TRUE)
#' }
#' @author Steven E. Pav \email{shabbychef@@gmail.com}
runExample <- function(port=NULL,launch.browser=TRUE,host=getOption('shiny.host','127.0.0.1'),display.mode='auto') { # nocov start
  appDir <- system.file("shiny-examples", "sadists", package = "sadists")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `sadists`.", call. = FALSE)
  }

  shiny::runApp(appDir, port=port, launch.browser=launch.browser,
								host=host, display.mode=display.mode)
} # nocov end

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
