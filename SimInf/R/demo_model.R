## SimInf, a framework for stochastic disease spread simulations
## Copyright (C) 2015  Pavol Bauer
## Copyright (C) 2015  Stefan Engblom
## Copyright (C) 2015  Stefan Widgren
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

##' Generate a model for demonstration
##'
##' @param nodes Number of nodes in the model. Default is 1.
##' @param days Number of days to model. Initializes \code{tspan} to
##' \code{\{0, 1, ..., days - 1\}}. Default is 1000 days.
##' @param model The name of the model. Default is 'SISe'.
##' @return A model
##' @include scheduled_events.R
##' @include SISe.R
##' @include SISe_sp.R
##' @include SISe3.R
##' @include SISe_sp.R
##' @export
##' @examples
##' ## Create a 'SISe3' demo model with 1 node and
##' ## initialize it to run over 1000 days.
##' model <- demo_model(nodes = 1, days = 1000, model = "SISe3")
##' result <- run(model)
##' plot(result)
demo_model <- function(nodes = 1,
                       days = 1000,
                       model = c("SISe", "SISe3", "SISe_sp"))
{
    ## Check 'nodes' argument
    if (!is.numeric(nodes))
        stop("'nodes' must be numeric.")
    if (!identical(length(nodes), 1L))
        stop("Length of 'nodes' must be one.")
    if (!is_wholenumber(nodes))
        stop("'nodes' must be integer")
    nodes <- as.integer(nodes)
    if (nodes[1] < 1)
        stop("'nodes' must be >= 1")

    ## Check 'days' argument
    if (!is.numeric(days))
        stop("'days' must be numeric.")
    if (!identical(length(days), 1L))
        stop("Length of 'days' must be one.")
    if (!is_wholenumber(days))
        stop("'days' must be integer")
    days <- as.integer(days)
    if (days[1] < 1)
        stop("'days' must be >= 1")

    ## Check 'model' argument
    model <- match.arg(model)

    if (identical(model, "SISe")) {
        u0 <- data.frame(S = rep(99, nodes), I = rep(1, nodes))

        model <- SISe(u0      = u0,
                      tspan   = seq_len(days) - 1,
                      events  = NULL,
                      phi     = rep(1, nodes),
                      upsilon = 0.017,
                      gamma   = 0.1,
                      alpha   = 1,
                      beta_t1 = 0.19,
                      beta_t2 = 0.085,
                      beta_t3 = 0.075,
                      beta_t4 = 0.185,
                      end_t1  = 91,
                      end_t2  = 182,
                      end_t3  = 273,
                      end_t4  = 365,
                      epsilon = 0.000011)
    } else if (identical(model, "SISe3")) {
        u0 <- data.frame(S_1 = rep(10, nodes), I_1 = rep( 0, nodes),
                         S_2 = rep(20, nodes), I_2 = rep( 0, nodes),
                         S_3 = rep(70, nodes), I_3 = rep( 0, nodes))

        model <- SISe3(u0        = u0,
                       tspan     = seq_len(days) - 1,
                       events    = NULL,
                       phi       = rep(1, nodes),
                       upsilon_1 = 0.0357,
                       upsilon_2 = 0.0357,
                       upsilon_3 = 0.00935,
                       gamma_1   = 0.1,
                       gamma_2   = 0.1,
                       gamma_3   = 0.1,
                       alpha     = 1.0,
                       beta_t1   = 0.19,
                       beta_t2   = 0.085,
                       beta_t3   = 0.075,
                       beta_t4   = 0.185,
                       end_t1    = 91,
                       end_t2    = 182,
                       end_t3    = 273,
                       end_t4    = 365,
                       epsilon   = 0.000011)
    } else if (identical(model, "SISe_sp")) {
        ## Place nodes on a grid
        nodes <- sqrt(nodes)
        if (!is_wholenumber(nodes))
            stop("'sqrt(nodes)' must be integer")
        distance <- expand.grid(x = seq_len(nodes),
                                y = seq_len(nodes))
        distance <- distance_matrix(distance$x, distance$y, 2)

        u0 <- data.frame(S = rep(99, nrow(distance)),
                         I = rep( 1, nrow(distance)))

        model <- SISe_sp(u0      = u0,
                         tspan   = seq_len(days) - 1,
                         events  = NULL,
                         phi     = rep(1, nrow(distance)),
                         upsilon = 0.017,
                         gamma   = 0.1,
                         alpha   = 1,
                         beta_t1 = 0.19,
                         beta_t2 = 0.085,
                         beta_t3 = 0.075,
                         beta_t4 = 0.185,
                         end_t1  = 91,
                         end_t2  = 182,
                         end_t3  = 273,
                         end_t4  = 365,
                         epsilon = 0.000011,
                         distance  = distance,
                         coupling  = 0.002)
    }

    model
}
