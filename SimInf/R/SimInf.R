## SimInf, a framework for stochastic disease spread simulations
## Copyright (C) 2015  Pavol Bauer
## Copyright (C) 2015 - 2016  Stefan Engblom
## Copyright (C) 2015 - 2016  Stefan Widgren
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.

##' A Framework for Stochastic Disease Spread Simulations
##'
##' @docType package
##' @name SimInf
##' @import methods
##' @useDynLib SimInf, .registration=TRUE
NULL

##' Unload hook function
##'
##' @param libpath A character string giving the complete path to the
##' package.
##' @keywords internal
.onUnload <- function (libpath)
{
    library.dynam.unload("SimInf", libpath)
}

##' Scheduled events example data
##'
##' Synthetic scheduled events data to demonstrate the \code{SISe3}
##' model. The data contains 783773 events for 1600 nodes over 365 * 4
##' days.
##' @name events_SISe3
##' @docType data
##' @usage data(events_SISe3)
##' @format A \code{data.frame}
##' @keywords dataset
##' @examples
##' \dontrun{
##' data(u0_SISe3)
##' data(events_SISe3)
##'
##' ## Initialize a model with 1600 nodes. Randomly select nodes to be
##' ## infected with a probability equal to 0.1.
##' model <- SISe3(u0        = u0_SISe3,
##'                tspan     = seq(0, 1460, by = 7),
##'                phi       = sample(0:1, 1600, TRUE, c(0.9, 0.1)),
##'                events    = events_SISe3,
##'                upsilon_1 = 0.04,
##'                upsilon_2 = 0.04,
##'                upsilon_3 = 0.01,
##'                gamma_1   = 0.1,
##'                gamma_2   = 0.1,
##'                gamma_3   = 0.1,
##'                alpha     = 1,
##'                beta_t1   = 0.15,
##'                beta_t2   = 0.10,
##'                beta_t3   = 0.10,
##'                beta_t4   = 0.15,
##'                end_t1    = 91,
##'                end_t2    = 182,
##'                end_t3    = 273,
##'                end_t4    = 365,
##'                epsilon   = 0)
##'
##' plot(model@events)
##'
##' result <- run(model)
##' plot(result)
##' }
NULL

##' Example data to initialize a model
##'
##' Synthetic init data for 1600 nodes to demonstrate the \code{SISe3}
##' model.
##' @name u0_SISe3
##' @docType data
##' @usage data(u0_SISe3)
##' @format A \code{data.frame}
##' @keywords dataset
##' @examples
##' \dontrun{
##' data(u0_SISe3)
##' data(events_SISe3)
##'
##' ## Initialize a model with 1600 nodes. Randomly select nodes to be
##' ## infected with a probability equal to 0.1.
##' model <- SISe3(u0        = u0_SISe3,
##'                tspan     = seq(0, 1460, by = 7),
##'                phi       = sample(0:1, 1600, TRUE, c(0.9, 0.1)),
##'                events    = events_SISe3,
##'                upsilon_1 = 0.04,
##'                upsilon_2 = 0.04,
##'                upsilon_3 = 0.01,
##'                gamma_1   = 0.1,
##'                gamma_2   = 0.1,
##'                gamma_3   = 0.1,
##'                alpha     = 1,
##'                beta_t1   = 0.15,
##'                beta_t2   = 0.10,
##'                beta_t3   = 0.10,
##'                beta_t4   = 0.15,
##'                end_t1    = 91,
##'                end_t2    = 182,
##'                end_t3    = 273,
##'                end_t4    = 365,
##'                epsilon   = 0)
##'
##' plot(model@events)
##'
##' result <- run(model)
##' plot(result)
##' }
NULL
