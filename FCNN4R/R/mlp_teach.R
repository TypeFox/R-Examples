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


#' Backpropagation (batch) teaching
#'
#' Backpropagation (a teaching algorithm) is a simple steepest
#' descent algorithm for MSE minimisation, in which weights are updated according
#' to (scaled) gradient of MSE.
#'
#' @note The name `backpropagation' is commonly used in two contexts, which
#' sometimes causes confusion. Firstly, backpropagation can be understood as
#' an efficient algorithm for MSE gradient computation that was first described
#' by Bryson and Ho in the '60s of 20th century and reinvented in the '80s.
#' Secondly, the name backpropagation is (more often) used to refer to the steepest
#' descent method that uses gradient of MSE computed efficiently by means
#' of the aforementioned algorithm. This ambiguity is probably caused by the fact
#' that in practically all neural network implementations, the derivatives of MSE
#' and weight updates are computed simultaneously in one backward pass (from
#' output layer to input layer).
#'
#' @param net an object of \code{mlp_net} class
#' @param input numeric matrix, each row corresponds to one input vector,
#'        the number of columns must be equal to the number of neurons
#'        in the network input layer
#' @param output numeric matrix with rows corresponding to expected outputs,
#'        the number of columns must be equal to the number of neurons
#'        in the network output layer, the number of rows must be equal to the number
#'        of input rows
#' @param tol_level numeric value, error (MSE) tolerance level
#' @param max_epochs integer value, maximal number of epochs (iterations)
#' @param learn_rate numeric value, learning rate in the backpropagation
#'        algorithm (default 0.7)
#' @param l2reg numeric value, L2 regularization parameter (default 0)
#' @param report_freq integer value, progress report frequency, if set to 0 no information is printed
#'        on the console (this is the default)
#'
#' @return Two-element list, the first field (\code{net}) contains the trained network,
#'         the second (\code{mse}) - the learning history (MSE in consecutive epochs).
#'
#' @references
#' A.E. Bryson and Y.C. Ho. \emph{Applied optimal control: optimization, estimation,
#' and control. Blaisdell book in the pure and applied sciences.} Blaisdell Pub. Co., 1969.
#'
#' David E. Rumelhart, Geoffrey E. Hinton, and Ronald J. Williams. \emph{Learning representations
#' by back-propagating errors.} Nature, 323(6088):533-536, October 1986.
#'
#'
#' @keywords teaching
#'
#' @export mlp_teach_bp
#'
mlp_teach_bp <- function(net, input, output,
                         tol_level, max_epochs,
                         learn_rate = 0.7, l2reg = 0,
                         report_freq = 0)
{
    if (tol_level <= 0) stop("tolerance level should be positive")
    if (learn_rate <= 0) stop("learning rate should be positive")
    if (l2reg < 0) stop("L2 regularization parameter should be nonnegative")

    gm <- mlp_grad(net, input, output)
    g <- gm$grad
    mse <- gm$mse
    if (mse < tol_level) {
        return(list(net = net, mse = NULL))
    }
    w0 <- mlp_get_weights(net)
    g <- g + l2reg * w0

    mseh <- numeric(length = max_epochs)
    for (i in 1:max_epochs) {
        w1 <- w0 - learn_rate * g
        net <- mlp_set_weights(net, w1)
        gm <- mlp_grad(net, input, output)
        g <- gm$grad + l2reg * w1
        mse <- gm$mse
        mseh[i] <- mse
        if (report_freq) {
            if (!(i %% report_freq)) {
                mes <- paste0("backpropagation; epoch ", i,
                              ", mse: ", mse, " (desired: ", tol_level, ")\n")
                cat(mes);
            }
        }
        if (mse < tol_level) break;
        w0 <- w1
    }
    if (mse > tol_level) {
        warning(paste0("algorithm did not converge, mse after ", i,
                       " epochs is ", mse, " (desired: ", tol_level, ")"))
    }
    return(list(net = net, mse = mseh[1:i]))
}




#' Rprop teaching
#'
#' Rprop is a fast and robust adaptive step method based on backpropagation.
#' For details, please refer to the original paper given in References section.
#'
#'
#' @param net an object of \code{mlp_net} class
#' @param input numeric matrix, each row corresponds to one input vector,
#'        the number of columns must be equal to the number of neurons
#'        in the network input layer
#' @param output numeric matrix with rows corresponding to expected outputs,
#'        the number of columns must be equal to the number of neurons
#'        in the network output layer, the number of rows must be equal to the number
#'        of input rows
#' @param tol_level numeric value, error (MSE) tolerance level
#' @param max_epochs integer value, maximal number of epochs (iterations)
#' @param l2reg numeric value, L2 regularization parameter (default 0)
#' @param u numeric value, Rprop algorithm parameter (default 1.2)
#' @param d numeric value, Rprop algorithm parameter (default 0.5)
#' @param gmax numeric value, Rprop algorithm parameter (default 50)
#' @param gmin numeric value, Rprop algorithm parameter (default 1e-6)
#' @param report_freq integer value, progress report frequency, if set to 0 no information is printed
#'        on the console (this is the default)
#'
#' @return Two-element list, the first field (\code{net}) contains the trained network,
#'         the second (\code{mse}) - the learning history (MSE in consecutive epochs).
#'
#' @references
#' M. Riedmiller. \emph{Rprop - Description and Implementation Details: Technical Report.} Inst. f.
#' Logik, Komplexitat u. Deduktionssysteme, 1994.
#'
#' @keywords teaching
#'
#' @export mlp_teach_rprop
#'
mlp_teach_rprop <- function(net, input, output,
                            tol_level, max_epochs,
                            l2reg = 0, u = 1.2, d = 0.5, gmax = 50., gmin = 1e-6,
                            report_freq = 0)
{
    if (tol_level <= 0) stop("tolerance level should be positive")
    if (l2reg < 0) stop("L2 regularization parameter should be nonnegative")

    gm <- mlp_grad(net, input, output)
    w0 <- mlp_get_weights(net);
    g0 <- gm$grad + l2reg * w0
    mse <- gm$mse
    if (mse < tol_level) {
        return(list(net = net, mse = NULL))
    }
    w1 <- w0 - 0.7 * g0;
    net <- mlp_set_weights(net, w1)

    mseh <- numeric(length = max_epochs)
    gm <- mlp_grad(net, input, output)
    g1 <- gm$grad + l2reg * w1
    mse <- gm$mse
    mseh[1] <- mse
    if (report_freq == 1) {
        mes <- paste0("Rprop; epoch 1",
                        ", mse: ", mse, " (desired: ", tol_level, ")\n")
        cat(mes);
    }
    if (mse < tol_level) {
        return(list(net = net, mse = mse))
    }

    nw <- length(w0)
    if (gmin > 1e-1) {
        gam <- gmin
    } else {
        gam <- min(0.1, gmax)
    }
    gamma <- rep(gam, nw)

    for (i in 2:max_epochs) {
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
        gm <- mlp_grad(net, input, output)
        g1 <- gm$grad + l2reg * w1
        mse <- gm$mse
        mseh[i] <- mse
        if (report_freq) {
            if (!(i %% report_freq)) {
                mes <- paste0("Rprop; epoch ", i, ", mse: ", mse,
                              " (desired: ", tol_level, ")\n")
                cat(mes);
            }
        }
        if (mse < tol_level) break;
    }
    if (mse > tol_level) {
        warning(paste0("algorithm did not converge, mse after ", i,
                       " epochs is ", mse, " (desired: ", tol_level, ")"))
    }
    return(list(net = net, mse = mseh[1:i]))
}



#' Stochastic gradient descent with (optional) RMS weights scaling, weight
#' decay, and momentum
#'
#' This function implements the stochastic gradient descent method with
#' optional modifications: L2 regularization, root mean square gradient scaling, weight decay,
#' and momentum.
#'
#' @param net an object of \code{mlp_net} class
#' @param input numeric matrix, each row corresponds to one input vector
#'        number of columns must be equal to the number of neurons
#'        in the network input layer
#' @param output numeric matrix with rows corresponding to expected outputs,
#'        number of columns must be equal to the number of neurons
#'        in the network output layer, number of rows must be equal to the number
#'        of input rows
#' @param tol_level numeric value, error (MSE) tolerance level
#' @param max_epochs integer value, maximal number of epochs (iterations)
#' @param learn_rate numeric value, (initial) learning rate, depending
#'                   on the problem at hand, learning rates of 0.001 or 0.01 should
#'                   give satisfactory convergence
#' @param l2reg numeric value, L2 regularization parameter (default 0)
#' @param minibatchsz integer value, the size of the mini batch (default 100)
#' @param lambda numeric value, rmsprop parameter controlling the update
#'               of mean squared gradient, reasonable value is 0.1 (default 0)
#' @param gamma numeric value, weight decay parameter (default 0)
#' @param momentum numeric value, momentum parameter, reasonable values are
#'                 between 0.5 and 0.9 (default 0)
#' @param report_freq integer value, progress report frequency, if set to 0
#'        no information is printed on the console (this is the default)
#'
#' @return Two-element list, the first field (\code{net}) contains the trained network,
#'         the second (\code{mse}) - the learning history (MSE in consecutive epochs).
#'
#' @keywords teaching
#'
#' @export mlp_teach_sgd
#'
mlp_teach_sgd <- function(net, input, output, tol_level, max_epochs,
                          learn_rate, l2reg = 0,
                          minibatchsz = 100, lambda = 0, gamma = 0,
                          momentum = 0,
                          report_freq = 0)
{
    if (tol_level <= 0) stop("tolerance level should be positive")
    if (learn_rate <= 0) stop("learning rate should be positive")
    if (l2reg < 0) stop("L2 regularization parameter should be nonnegative")
    N <- dim(input)[1]
    M <- round(minibatchsz)
    if ((M < 1) || (M >= N)) {
        stop("minibatch size should be at least 1 and less than the number of records")
    }
    W <- mlp_get_no_active_w(net)
    if (lambda != 0) {
        ms <- rep(1, W)
    }
    if (momentum != 0) {
        mm <- rep(0, W)
    }
    idx <- sample.int(N, M)
    gm <- mlp_grad(net, input[idx, , drop = FALSE], output[idx, , drop = FALSE])
    w0 <- mlp_get_weights(net)
    g <- gm$grad + l2reg * w0
    mse <- gm$mse
    if (mse < tol_level) {
        if (mlp_mse(net, input, output) < tol_level) {
            return(list(net = net, mse = NULL))
        }
    }

    mseh <- numeric(length = max_epochs)
    for (i in 1:max_epochs) {
        dw <- -learn_rate * g
        if (lambda != 0) {
            dw <- dw / sqrt(ms)
            ms <- (1 - lambda) * ms + lambda * g^2
        }
        if (gamma != 0) {
            dw <- dw / (1 + gamma * (i - 1))
        }
        if (momentum != 0) {
            dw <- momentum * mm + dw
            mm <- dw
        }
        w1 <- w0 + dw
        net <- mlp_set_weights(net, w1)
        idx <- sample.int(N, M)
        gm <- mlp_grad(net, input[idx, , drop = FALSE], output[idx, , drop = FALSE])
        g <- gm$grad + l2reg * w1
        mse <- gm$mse
        mseall <- FALSE
        mseh[i] <- mse
        if (report_freq) {
            if (!(i %% report_freq)) {
                mes <- paste0("stochastic gradient descent; epoch ", i,
                              ", mse: ", mse, " (desired: ", tol_level, ")\n")
                cat(mes)
            }
        }
        # compute mse on the entire training set if mse on the minibatch
        # is less than tolerance level
        if (mse < tol_level) {
            mse <- mlp_mse(net, input, output)
            mseh[i] <- mse
            mseall <- TRUE
            if (mse < tol_level) break;
        }
        w0 <- w1
    }
    # compute mse on the entire training set
    if (!mseall) {
        mse <- mlp_mse(net, input, output)
        mseh[i] <- mse
    }
    if (mse > tol_level) {
        warning(paste0("algorithm did not converge, mse (on the entire training set) after ", i,
                       " epochs is ", mse, " (desired: ", tol_level, ")"))
    }
    return(list(net = net, mse = mseh[1:i]))
}
