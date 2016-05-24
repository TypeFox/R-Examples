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

##' Class \code{"SISe3"}
##'
##' Class to handle the SISe3 \code{\link{siminf_model}} model.
##' @include siminf_model.R
##' @include AllGenerics.R
##' @export
setClass("SISe3", contains = c("siminf_model"))

##' Create a SISe3 model
##'
##' Create a SISe3 model to be used by the simulation framework.
##'
##'
##' The argument \code{u0} must be a \code{data.frame} with one row for
##' each node with the following columns:
##' \describe{
##' \item{S_1}{The number of sucsceptible in age category 1}
##' \item{I_1}{The number of infected in age category 1}
##' \item{S_2}{The number of sucsceptible in age category 2}
##' \item{I_2}{The number of infected in age category 2}
##' \item{S_3}{The number of sucsceptible in age category 3}
##' \item{I_3}{The number of infected in age category 3}
##' }
##'
##' @template beta-section
##' @param u0 A \code{data.frame} with the initial state in each
##' node, see details.
##' @param tspan An increasing sequence of points in time where the
##' state of the system is to be returned.
##' @param events a \code{data.frame} with the scheduled events, see
##' \code{\link{siminf_model}}.
##' @param phi A numeric vector with the initial environmental
##' infectious pressure in each node. Default NULL which gives 0 in
##' each node.
##' @param upsilon_1 Indirect transmission rate of the environmental
##' infectious pressure in age category 1
##' @param upsilon_2 Indirect transmission rate of the environmental
##' infectious pressure in age category 2
##' @param upsilon_3 Indirect transmission rate of the environmental
##' infectious pressure in age category 3
##' @param gamma_1 The recovery rate from infected to susceptible for
##' age category 1
##' @param gamma_2 The recovery rate from infected to susceptible for
##' age category 2
##' @param gamma_3 The recovery rate from infected to susceptible for
##' age category 3
##' @param alpha Shed rate from infected individuals
##' @template beta-param
##' @param epsilon The background environmental infectious pressure
##' @return \code{SISe3}
##' @include check_arguments.R
##' @export
SISe3 <- function(u0,
                  tspan,
                  events    = NULL,
                  phi       = NULL,
                  upsilon_1 = NULL,
                  upsilon_2 = NULL,
                  upsilon_3 = NULL,
                  gamma_1   = NULL,
                  gamma_2   = NULL,
                  gamma_3   = NULL,
                  alpha     = NULL,
                  beta_t1   = NULL,
                  beta_t2   = NULL,
                  beta_t3   = NULL,
                  beta_t4   = NULL,
                  end_t1    = NULL,
                  end_t2    = NULL,
                  end_t3    = NULL,
                  end_t4    = NULL,
                  epsilon   = NULL)
{
    compartments <- c("S_1", "I_1", "S_2", "I_2", "S_3", "I_3")

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

    ## Check 'gdata' parameters
    check_gdata_arg(upsilon_1, upsilon_2, upsilon_3, gamma_1, gamma_2, gamma_3,
                    alpha, beta_t1, beta_t2, beta_t3, beta_t4, epsilon)

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

    E <- Matrix(c(1, 0, 0, 1, 0, 0,
                  0, 0, 0, 1, 0, 0,
                  0, 1, 0, 0, 1, 0,
                  0, 0, 0, 0, 1, 0,
                  0, 0, 1, 0, 0, 1,
                  0, 0, 0, 0, 0, 1),
                nrow   = 6,
                ncol   = 6,
                byrow  = TRUE,
                sparse = TRUE)
    E <- as(E, "dgCMatrix")
    colnames(E) <- as.character(1:6)
    rownames(E) <- compartments

    N <- matrix(c(2, 0,
                  2, 0,
                  0, 2,
                  0, 2,
                  0, 0,
                  0, 0),
                nrow   = 6,
                ncol   = 2,
                byrow  = TRUE)
    colnames(N) <- as.character(1:2)
    rownames(N) <- compartments

    G <- Matrix(c(1, 1, 0, 0, 0, 0,
                  1, 1, 0, 0, 0, 0,
                  0, 0, 1, 1, 0, 0,
                  0, 0, 1, 1, 0, 0,
                  0, 0, 0, 0, 1, 1,
                  0, 0, 0, 0, 1, 1),
                nrow   = 6,
                ncol   = 6,
                byrow  = TRUE,
                sparse = TRUE)
    G <- as(G, "dgCMatrix")
    colnames(G) <- as.character(1:6)
    rownames(G) <- c("S_1 -> I_1", "I_1 -> S_1",
                     "S_2 -> I_2", "I_2 -> S_2",
                     "S_3 -> I_3", "I_3 -> S_3")

    S <- Matrix(c(-1,  1,  0,  0,  0,  0,
                   1, -1,  0,  0,  0,  0,
                   0,  0, -1,  1,  0,  0,
                   0,  0,  1, -1,  0,  0,
                   0,  0,  0,  0, -1,  1,
                   0,  0,  0,  0,  1, -1),
                nrow   = 6,
                ncol   = 6,
                byrow  = TRUE,
                sparse = TRUE)
    S <- as(S, "dgCMatrix")
    colnames(S) <- as.character(1:6)
    rownames(S) <- compartments

    v0 <- matrix(phi, nrow  = 1, byrow = TRUE)
    storage.mode(v0) <- "double"

    ldata <- matrix(c(end_t1, end_t2, end_t3, end_t4),
                    nrow  = 4,
                    byrow = TRUE)
    storage.mode(ldata) <- "double"

    gdata <- c(upsilon_1, upsilon_2, upsilon_3,
               gamma_1, gamma_2, gamma_3,
               alpha,
               beta_t1, beta_t2, beta_t3, beta_t4,
               epsilon)
    storage.mode(gdata) <- "double"
    names(gdata) <- c("upsilon_1", "upsilon_2", "upsilon_3",
                      "gamma_1", "gamma_2", "gamma_3",
                      "alpha",
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

    return(as(model, "SISe3"))
}

##' @rdname susceptible-methods
##' @export
setMethod("susceptible",
          signature("SISe3"),
          function(model, age = 1:3, i = NULL, by = 1, ...)
          {
              if (identical(dim(model@U), c(0L, 0L)))
                  stop("Please run the model first, the 'U' matrix is empty")

              age_categories <- 1:3
              stopifnot(all(age %in% age_categories))

              result <- NULL
              j <- seq(from = 1, to = dim(model@U)[2], by = by)

              for (k in age_categories) {
                  ## Are we interested in this age category?
                  if (k %in% age) {
                      ## Select rows for the specific age category
                      ii <- seq(from = 1 + (k - 1) * 2, to = dim(model@U)[1], by = 6)

                      ## Are we only interested in susceptible from
                      ## specific nodes?
                      if (!is.null(i))
                          ii <- ii[i]

                      ## Extract susceptible and add to result
                      if (is.null(result)) {
                          result <- as.matrix(model@U[ii, j, drop = FALSE])
                      } else {
                          result <- result + as.matrix(model@U[ii, j, drop = FALSE])
                      }
                  }
              }

              result
          }
)

##' @rdname infected-methods
##' @export
setMethod("infected",
          signature("SISe3"),
          function(model, age = 1:3, i = NULL, by = 1, ...)
          {
              if (identical(dim(model@U), c(0L, 0L)))
                  stop("Please run the model first, the 'U' matrix is empty")

              age_categories <- 1:3
              stopifnot(all(age %in% age_categories))

              result <- NULL
              j <- seq(from = 1, to = dim(model@U)[2], by = by)

              for (k in age_categories) {
                  ## Are we interested in this age category?
                  if (k %in% age) {
                      ## Select rows for the specific age category
                      ii <- seq(from = k * 2, to = dim(model@U)[1], by = 6)

                      ## Are we only interested in infected from
                      ## specific nodes?
                      if (!is.null(i))
                          ii <- ii[i]

                      ## Extract infected and add to result
                      if (is.null(result)) {
                          result <- as.matrix(model@U[ii, j, drop = FALSE])
                      } else {
                          result <- result + as.matrix(model@U[ii, j, drop = FALSE])
                      }
                  }
              }

              result
          }
)

##' @rdname prevalence-methods
##' @export
setMethod("prevalence",
          signature("SISe3"),
          function(model, age = 1:3, wnp = FALSE, i = NULL, by = 1, ...)
          {
              I <- infected(model, age = age, i = i, by = by)
              S <- susceptible(model, age = age, i = i, by = by)

              if (identical(wnp, FALSE)) {
                  I <- colSums(I)
                  S <- colSums(S)
              }

              I / (S + I)
          }
)

##' @name plot-methods
##' @aliases plot plot-methods plot,SISe3-method
##' @docType methods
##' @importFrom graphics plot
##' @export
setMethod("plot",
          signature(x = "SISe3"),
          function(x, t0 = NULL, ...)
      {
          callNextMethod(x,
                         t0 = t0,
                         legend = expression(S[1], I[1], S[2], I[2], S[3], I[3]),
                         ...)
      }
)
