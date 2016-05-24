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


#' Minimum magnitude pruning
#'
#' Minimum magnitude pruning is a brute force, easy-to-implement pruning algorithm
#' in which in each step the weight with the smallest absolute value is turned off. This
#' algorithm requires reteaching network in almost every step and yields suboptimal results.
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
#' @param max_reteach_epochs integer value, maximal number of epochs (iterations) allowed
#'        when reteaching network
#' @param report logical value, if TRUE, information about the pruning process
#'        will be printed on the console (FALSE by default)
#' @param plots logical value, if TRUE, the initial network is plotted and then
#'        replotted every time neuron is removed and at the end of pruning (FALSE by default)
#'
#' @return Three-element list, the first field (\code{net}) contains the pruned network,
#'         the second (\code{wcount}) - the number of connections removed (inactivated),
#'         the third (\code{ncount}) - the number of neurons removed.
#'
#' @keywords pruning
#'
#' @export mlp_prune_mag
#'
mlp_prune_mag <- function(net, input, output,
                          tol_level, max_reteach_epochs,
                          report, plots = FALSE)
{
    if (tol_level <= 0.) {
        stop("tolerance level should be positive");
    }
    mse <- mlp_mse(net, input, output)
    if (mse > tol_level) {
        stop(paste0("network should be trained with MSE reduced to given tolerance level (",
                    tol_level, ") before pruning; MSE is ", mse))
    }

    stop <- FALSE;
    countw <- 0
    countn <- 0

    if (plots) {
        mlp_plot(net)
    }

    while (!stop) {
        wts <- mlp_get_weights(net)
        wi <- which.min(abs(wts))
        wi <- mlp_get_w_abs_idx(net, wi)
        net <- mlp_set_w_st(net, FALSE, idx = wi)
        countw <- countw + 1;

        mse <- mlp_mse(net, input, output)
        if (mse > tol_level) {
            retres <- suppressWarnings(mlp_teach_rprop(net, input, output,
                                                       tol_level,
                                                       max_reteach_epochs))
            mse <- retres[[2]][length(retres[[2]])]
            if (mse > tol_level) {
                stop <- TRUE;
                countw <- countw - 1;
                net <- mlp_set_w_st(net, TRUE, idx = wi)
                net <- mlp_set_weights(net, wts)
            } else {
                net <- retres$net
                if (report) {
                    cat(paste0("removed weight ", wi, " from ", mlp_get_no_w(net),
                               " total (", mlp_get_no_active_w(net),
                               " remain active); network has been retrained\n"))
                }
            }
        } else {
            if (report) {
                cat(paste0("removed weight ", wi, " from ", mlp_get_no_w(net),
                            " total (", mlp_get_no_active_w(net),
                            " remain active)\n"))
            }
        }
        rmnres <- mlp_rm_neurons(net, report = report)
        net <- rmnres$net
        countn <- countn + rmnres$ncount
        countw <- countw + rmnres$wcount
        if ((rmnres$ncount != 0) && plots) {
            mlp_plot(net)
        }
    }
    if (report) {
        cat(paste0("pruning stopped, removed ", countw, " weight(s) and ", countn, " neuron(s)\n"))
    }
    if (plots) {
        mlp_plot(net)
    }

    return(list(net = net, wcount = countw, ncount = countn))
}


#' Optimal Brain Surgeon pruning
#'
#' The Optimal Brain Surgeon algorithm is a robust (yet computationally demanding)
#' pruning algorithm in which candidate weight to be turned off is determined
#' based on information about the inverse of (approximate) Hessian matrix of the MSE.
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
#' @param max_reteach_epochs integer value, maximal number of epochs (iterations) allowed
#'        when reteaching network
#' @param report logical value, if TRUE, information about the pruning process
#'        will be printed on the console (FALSE by default)
#' @param plots logical value, if TRUE, the initial network is plotted and then
#'        replotted every time neuron is removed and at the end of pruning (FALSE by default)
#' @param alpha numeric value, scaling factor used for initial Hessian approximation
#'
#' @return Three-element list, the first field (\code{net}) contains the pruned network,
#'         the second (\code{wcount}) - the number of connections removed (inactivated),
#'         the third (\code{ncount}) - the number of neurons removed.
#'
#' @references
#' B. Hassibi, D. G. Stork, and G. J. Wolff. \emph{Optimal Brain Surgeon
#' and General Network Pruning.} Technical Report CRC-TR-9235, RICOH California
#' Research Centre, 1992.
#'
#' @keywords pruning
#'
#' @export mlp_prune_obs
#'
mlp_prune_obs <- function(net, input, output,
                          tol_level, max_reteach_epochs,
                          report, plots = FALSE, alpha = 1e-5)
{
    if (tol_level <= 0.) {
        stop("tolerance level should be positive");
    }
    mse <- mlp_mse(net, input, output)
    if (mse > tol_level) {
        stop(paste0("network should be trained with MSE reduced to given tolerance level (",
                    tol_level, ") before pruning; MSE is ", mse))
    }

    stop <- FALSE;
    countw <- 0
    countn <- 0
    PN <- dim(output)
    P <- PN[1]
    N <- PN[2]
    PN <- prod(PN)

    if (plots) {
        mlp_plot(net)
    }

    while (!stop) {
        W <- mlp_get_no_active_w(net)
        H <- diag(1 / alpha, nrow = W, ncol = W)
        for (i in 1:P) {
            grads <- mlp_gradij(net, input, i)
            H <- .C("ihessupdate", as.integer(W), as.integer(N), as.numeric(PN),
                    grads, res = H)$res
        }
        wts <- mlp_get_weights(net)
        L <- .5 * wts^2 / diag(H)
        wi <- which.min(L)
        dw <- wts[wi] * H[, wi] / as.numeric(H[wi, wi])
        net <- mlp_set_weights(net, wts - dw)
        wi <- mlp_get_w_abs_idx(net, wi)
        net <- mlp_set_w_st(net, FALSE, idx = wi)
        countw <- countw + 1;

        mse <- mlp_mse(net, input, output)
        if (mse > tol_level) {
            retres <- suppressWarnings(mlp_teach_rprop(net, input, output,
                                                       tol_level,
                                                       max_reteach_epochs))
            mse <- retres[[2]][length(retres[[2]])]
            if (mse > tol_level) {
                stop <- TRUE;
                countw <- countw - 1;
                net <- mlp_set_w_st(net, TRUE, idx = wi)
                net <- mlp_set_weights(net, wts)
            } else {
                net <- retres$net
                if (report) {
                    cat(paste0("removed weight ", wi, " from ", mlp_get_no_w(net),
                               " total (", mlp_get_no_active_w(net),
                               " remain active); network has been retrained\n"))
                }
            }
        } else {
            if (report) {
                cat(paste0("removed weight ", wi, " from ", mlp_get_no_w(net),
                            " total (", mlp_get_no_active_w(net),
                            " remain active)\n"))
            }
        }
        rmnres <- mlp_rm_neurons(net, report = report)
        net <- rmnres$net
        countn <- countn + rmnres$ncount
        countw <- countw + rmnres$wcount
        if ((rmnres$ncount != 0) && plots) {
            mlp_plot(net)
        }
    }
    if (report) {
        cat(paste0("pruning stopped, removed ", countw, " weight(s) and ", countn, " neuron(s)\n"))
    }
    if (plots) {
        mlp_plot(net)
    }

    return(list(net = net, wcount = countw, ncount = countn))
}



