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

##' Class \code{"SISe"}
##'
##' Class to handle the SISe \code{\link{siminf_model}}.
##' @include siminf_model.R
##' @include AllGenerics.R
##' @export
setClass("SISe", contains = c("siminf_model"))

##' Create a SISe model
##'
##' Create a SISe model to be used by the simulation framework.
##'
##'
##' The argument \code{u0} must be a \code{data.frame} with one row for
##' each node with the following columns:
##' \describe{
##' \item{S}{The number of sucsceptible in each node}
##' \item{I}{The number of infected in each node}
##' }
##'
##' @template beta-section
##' @param u0 A \code{data.frame} with the initial state in each node,
##'     see details.
##' @param tspan An increasing sequence of points in time where the
##'     state of the system is to be returned.
##' @param events a \code{data.frame} with the scheduled events, see
##'     \code{\link{siminf_model}}.
##' @param phi A numeric vector with the initial environmental
##'     infectious pressure in each node. Default NULL which gives 0
##'     in each node.
##' @param upsilon Indirect transmission rate of the environmental
##'     infectious pressure
##' @param gamma The recovery rate from infected to susceptible
##' @param alpha Shed rate from infected individuals
##' @template beta-param
##' @param epsilon The background environmental infectious pressure
##' @return \code{SISe}
##' @include check_arguments.R
##' @export
SISe <- function(u0,
                 tspan,
                 events  = NULL,
                 phi     = NULL,
                 upsilon = NULL,
                 gamma   = NULL,
                 alpha   = NULL,
                 beta_t1 = NULL,
                 beta_t2 = NULL,
                 beta_t3 = NULL,
                 beta_t4 = NULL,
                 end_t1  = NULL,
                 end_t2  = NULL,
                 end_t3  = NULL,
                 end_t4  = NULL,
                 epsilon = NULL)
{
    compartments <- c("S", "I")

    ## Check arguments.

    ## Check u0
    if (!is.data.frame(u0))
        stop("'u0' must be a data.frame")
    if (!all(compartments %in% names(u0)))
        stop("Missing columns in u0")
    u0 <- u0[, compartments]

    ## Check initial infectious pressure
    if (is.null(phi))
        phi <- rep(0, nrow(u0))
    check_infectious_pressure_arg(nrow(u0), phi)

    ## Check for non-numeric parameters
    check_gdata_arg(upsilon, gamma, alpha, beta_t1, beta_t2, beta_t3, beta_t4,
                    epsilon)

    ## Check interval endpoints
    check_integer_arg(end_t1, end_t2, end_t3, end_t4)
    if (identical(length(end_t1), 1L))
        end_t1 <- rep(end_t1, nrow(u0))
    if (identical(length(end_t2), 1L))
        end_t2 <- rep(end_t2, nrow(u0))
    if (identical(length(end_t3), 1L))
        end_t3 <- rep(end_t3, nrow(u0))
    if (identical(length(end_t4), 1L))
        end_t4 <- rep(end_t4, nrow(u0))
    check_end_t_arg(nrow(u0), end_t1, end_t2, end_t3, end_t4)

    ## Arguments seems ok...go on

    E <- Matrix(c(1, 1,
                  0, 1),
                nrow   = 2,
                ncol   = 2,
                byrow  = TRUE,
                sparse = TRUE)
    E <- as(E, "dgCMatrix")
    colnames(E) <- as.character(1:2)
    rownames(E) <- compartments

    N <- matrix(integer(0), nrow = 0, ncol = 0)

    G <- Matrix(c(1, 1,
                  1, 1),
                nrow = 2,
                ncol = 2,
                byrow  = TRUE,
                sparse = TRUE)
    G <- as(G, "dgCMatrix")
    colnames(G) <- as.character(1:2)
    rownames(G) <- c("S -> I", "I -> S")

    S <- Matrix(c(-1,  1,
                   1, -1),
                nrow   = 2,
                ncol   = 2,
                byrow  = TRUE,
                sparse = TRUE)
    S <- as(S, "dgCMatrix")
    colnames(S) <- as.character(1:2)
    rownames(S) <- compartments

    v0 <- matrix(phi, nrow  = 1, byrow = TRUE)
    storage.mode(v0) <- "double"

    ldata <- matrix(c(end_t1, end_t2, end_t3, end_t4),
                    nrow  = 4,
                    byrow = TRUE)
    storage.mode(ldata) <- "double"

    gdata <- c(upsilon, gamma, alpha, beta_t1, beta_t2, beta_t3, beta_t4,
               epsilon)
    storage.mode(gdata) <- "double"
    names(gdata) <- c("upsilon", "gamma", "alpha",
                      "beta_t1", "beta_t2", "beta_t3", "beta_t4",
                      "epsilon")

    model <- siminf_model(G      = G,
                          S      = S,
                          E      = E,
                          N      = N,
                          tspan  = tspan,
                          events = events,
                          ldata  = ldata,
                          gdata  = gdata,
                          u0     = u0,
                          v0     = v0)

    return(as(model, "SISe"))
}

##' @rdname susceptible-methods
##' @export
setMethod("susceptible",
          signature("SISe"),
          function(model, i = NULL, by = 1, ...) {
              if (identical(dim(model@U), c(0L, 0L)))
                  stop("Please run the model first, the 'U' matrix is empty")

              ii <- seq(from = 1, to = dim(model@U)[1], by = 2)
              if (!is.null(i))
                  ii <- ii[i]
              j <- seq(from = 1, to = dim(model@U)[2], by = by)
              as.matrix(model@U[ii, j, drop = FALSE])
          }
)

##' @rdname infected-methods
##' @export
setMethod("infected",
          signature("SISe"),
          function(model, i = NULL, by = 1, ...) {
              if (identical(dim(model@U), c(0L, 0L)))
                  stop("Please run the model first, the 'U' matrix is empty")

              ii <- seq(from = 2, to = dim(model@U)[1], by = 2)
              if (!is.null(i))
                  ii <- ii[i]
              j <- seq(from = 1, to = dim(model@U)[2], by = by)
              as.matrix(model@U[ii, j, drop = FALSE])
          }
)

##' @rdname prevalence-methods
##' @export
setMethod("prevalence",
          signature("SISe"),
          function(model, wnp = FALSE, i = NULL, by = 1, ...) {
              I <- infected(model = model, i = i, by = by)
              S <- susceptible(model = model, i = i, by = by)

              if (identical(wnp, FALSE)) {
                  I <- colSums(I)
                  S <- colSums(S)
              }

              I / (S + I)
          }
)

##' @name plot-methods
##' @aliases plot plot-methods plot,SISe-method
##' @importFrom graphics plot
##' @export
setMethod("plot",
          signature(x = "SISe"),
          function(x, t0 = NULL, ...)
      {
          callNextMethod(x, t0 = t0, legend = c("S", "I"), ...)
      }
)
