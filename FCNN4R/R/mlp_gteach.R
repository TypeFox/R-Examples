# #########################################################################
# This file is a part of FCNN4R.
#
# Copyright (c) Grzegorz Klima 2015-2016
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
# #########################################################################


#' Rprop teaching - minimising arbitrary objective function
#'
#' This implementation (`generalisation') of the Rprop algorithm allows users to teach
#' network to minimise arbitrary objective function provided that functions
#' evaluating objective and computing gradient are provided.
#'
#' @param net an object of \code{mlp_net} class
#' @param obj_func function taking an object of \code{mlp_class} class
#'                 as a single argument returning objective to be minimised
#' @param gradient function taking an object of \code{mlp_class} class
#'                 as a single argument returning gradient of the objective
#' @param epochs integer value, number of epochs (iterations)
#' @param stop function (or NULL), a function taking objective history to date
#'             and returning Boolean value (if TRUE is returned, algorithm stops)
#'             (the default is not to stop until all iterations are performed)
#' @param report_freq integer value, progress report frequency, if set to 0
#'        no information is printed on the console (this is the default)
#' @param report_action function (or NULL), additional action to be taken while
#'                      printing progress reports, this should be a function
#'                      taking network as a single argument (default NULL)
#' @param u numeric value, Rprop algorithm parameter (default 1.2)
#' @param d numeric value, Rprop algorithm parameter (default 0.5)
#' @param gmax numeric value, Rprop algorithm parameter (default 50)
#' @param gmin numeric value, Rprop algorithm parameter (default 1e-6)
#'
#' @return Two-element list, the first field (\code{net}) contains the trained network,
#'         the second (\code{obj}) - the learning history (value of the objective
#'         function in consecutive epochs).
#'
#' @references
#' M. Riedmiller. \emph{Rprop - Description and Implementation Details: Technical Report.} Inst. f.
#' Logik, Komplexitat u. Deduktionssysteme, 1994.
#'
#' @examples
#' \dontrun{
#' # set up XOR problem
#' inp <- c(0, 0, 1, 1, 0, 1, 0, 1)
#' dim(inp) <- c(4, 2)
#' outp <- c(0, 1, 1, 0)
#' dim(outp) <- c(4, 1)
#' # objective
#' obj <- function(net)
#' {
#'     return(mlp_mse(net, inp, outp))
#' }
#' # gradient
#' grad <- function(net)
#' {
#'     return(mlp_grad(net, inp, outp)$grad)
#' }
#' # stopping citerion
#' tol <- function(oh) {
#'     if (oh[length(oh)] <= 5e-5) { return(TRUE); }
#'     return(FALSE)
#' }
#' # create a 2-6-1 network
#' net <- mlp_net(c(2, 6, 1))
#' # set activation function in all layers
#' net <- mlp_set_activation(net, layer = "a", "sigmoid")
#' # randomise weights
#' net <- mlp_rnd_weights(net)
#' # teach
#' netobj <- mlp_teach_grprop(net, obj, grad, epochs = 500,
#'                            stop = tol,
#'                            report_freq = 1)
#' # plot learning history
#' plot(netobj$obj, type = 'l')
#' }
#'
#' @keywords teaching
#'
#' @export mlp_teach_grprop
#'
mlp_teach_grprop <- function(net, obj_func, gradient,
                             epochs, stop = NULL,
                             report_freq = 0, report_action = NULL,
                             u = 1.2, d = 0.5, gmax = 50., gmin = 1e-6)
{
    if (!is.function(obj_func)) stop("obj_func should be a function")
    if (!is.function(gradient)) stop("gradient should be a function")
    if (is.null(stop)) {
        stop <- function(objh) { return(FALSE); }
    } else {
        if (!is.function(stop)) stop("stop should be a function (or NULL)")
    }

    g0 <- gradient(net)
    obj <- obj_func(net)
    w0 <- mlp_get_weights(net);
    w1 <- w0 - 0.7 * g0;
    net <- mlp_set_weights(net, w1)

    objh <- numeric(length = epochs)
    g1 <- gradient(net)
    obj <- obj_func(net)
    objh[1] <- obj
    if (report_freq == 1) {
        mes <- paste0("Rprop; epoch 1, objective: ", obj, "\n")
        cat(mes);
        if (!is.null(report_action)) {
            report_action(net)
        }
    }
    if (stop(objh[1:1])) return(list(net = net, obj = objh[1:1]));

    nw <- length(w0)
    if (gmin > 1e-1) {
        gam <- gmin
    } else {
        gam <- min(0.1, gmax)
    }
    gamma <- rep(gam, nw)

    for (i in 2:epochs) {
        # determine step and update gamma
        dw <- rep(0, nw)
        ig0 <- (g1 > 0)
        il0 <- (g1 < 0)
        i1 <- (g0 * g1 > 0)
        ind <- which(i1 & ig0)
        dw[ind] <- -gamma[ind]
        ind <- which(i1 & !ig0)
        dw[ind] <- gamma[ind]
        gamma[i1] <- pmin(u * gamma[i1], gmax)
        i2 <- (g0 * g1 < 0)
        ind <- which(i2)
        dw[ind] <- 0
        gamma[ind] <- pmax(d * gamma[ind], gmin)
        i3 <- (g0 * g1 == 0)
        ind <- which(i3 & ig0)
        dw[ind] <- -gamma[ind]
        ind <- which(i3 & il0)
        dw[ind] <- gamma[ind]
        # update weights
        w0 <- mlp_get_weights(net);
        w1 <- w0 + dw;
        net <- mlp_set_weights(net, w1)
        # new gradients
        g0 <- g1;
        g1 <- gradient(net)
        obj <- obj_func(net)
        objh[i] <- obj
        if (report_freq) {
            if (!(i %% report_freq)) {
                mes <- paste0("Rprop; epoch ", i, ", objective: ", obj, "\n")
                cat(mes);
                if (!is.null(report_action)) {
                    report_action(net)
                }
            }
        }
        if (stop(objh[1:i])) break;
    }
    return(list(net = net, obj = objh[1:i]))
}



#' Teaching networks using Simulated Annealing
#'
#' This function can be used to teach an ANN to minimise arbitrary objective
#' function.
#'
#' @param net an object of \code{mlp_net} class
#' @param obj_func function taking an object of \code{mlp_class} class
#'                 as a single argument returning objective to be minimised
#' @param Tinit numeric value, initial temperature (default is 1)
#' @param epochs integer value, number of epochs (iterations) (default is 1000)
#' @param report_freq integer value, progress report frequency, if set to 0
#'        no information is printed on the console (this is the default)
#' @param report_action function (or NULL), additional action to be taken while
#'                      printing progress reports, this should be a function
#'                      taking network as a single argument (default NULL)
#'
#' @return Two-element list, the first field (\code{net}) contains the trained network,
#'         the second (\code{obj}) - the learning history (value of the objective
#'         function in consecutive epochs).
#'
#' @examples
#' \dontrun{
#' # set up XOR problem
#' inp <- c(0, 0, 1, 1, 0, 1, 0, 1)
#' dim(inp) <- c(4, 2)
#' outp <- c(0, 1, 1, 0)
#' dim(outp) <- c(4, 1)
#' # objective
#' obj <- function(net)
#' {
#'     return(mlp_mse(net, inp, outp))
#' }
#' # create a 2-6-1 network
#' net <- mlp_net(c(2, 6, 1))
#' # set activation function in all layers
#' net <- mlp_set_activation(net, layer = "a", "sigmoid")
#' # teach
#' netobj <- mlp_teach_sa(net, obj, Tinit = 1, epochs = 1000,
#'                        report_freq = 1)
#' # plot learning history
#' plot(netobj$obj, type = 'l')
#' }
#'
#' @keywords teaching
#'
#' @export mlp_teach_sa
#'
mlp_teach_sa <- function(net, obj_func, Tinit = 1, epochs = 1000,
                         report_freq = 0, report_action = NULL)
{
    if (!is.function(obj_func)) stop("obj_func should be a function")
    if (Tinit <= 0) stop("initial temperature should be positive")
    obj <- obj_func(net)
    bf <- obj
    w <- mlp_get_weights(net)
    W <- length(w)
    best <- w
    objh <- numeric(length = epochs)
    fevalcount <- 1

    for (i in 1:epochs) {
        T <- Tinit * abs(i / epochs - 1) ^ 2
        repeat {
            # step
            dw <- numeric(W)
            noind <- 1 + floor(W * (1 - i / epochs))
            dw[sample(1:W, noind)] <- rnorm(noind) * .1
            w1 <- w + dw
            net <- mlp_set_weights(net, w1)
            obj1 <- obj_func(net)
            fevalcount <- fevalcount + 1
            # improvement?
            if (obj1 < obj) { break; }
            # acceptance probability
            pr <- exp(-1000 * (obj1 - obj) / T)
            # handle NaNs and Infs
            if (!is.finite(pr)) { pr <- 0; }
            # accept?
            if (runif(1) < pr) { break; }
        }
        w <- w1
        obj <- obj1
        objh[i] <- obj
        if (bf > obj) {
            bf <- obj
            best <- w
        }

        if (report_freq) {
            if (!(i %% report_freq)) {
                mes <- paste0("simulated annealing; epoch ", i,
                              " (obj. f. evaluations: ", fevalcount, ")\n",
                              "objective: ", obj)
                if (bf < obj) {
                    mes <- paste0(mes, " (best ever: ", bf, ")\n")
                } else {
                    mes <- paste0(mes, "\n")
                }
                cat(mes);
                if (!is.null(report_action)) {
                    report_action(net)
                }
            }
        }
    }
    if (report_freq) {
        mes <- paste0("simulated annealing stopped after ", epochs, " epochs (",
                    fevalcount, " objective evaluations),\nobjective: ", bf, "\n")
        cat(mes);
    }
    net <- mlp_set_weights(net, best)
    return(list(net = net, obj = objh[1:i]))
}



