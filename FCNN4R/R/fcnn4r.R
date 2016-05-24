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

.onAttach <- function(...)
{
  packageStartupMessage(paste0("Fast Compressed Neural Networks for R ",
                        utils::packageVersion("FCNN4R"), " (FCNN library ver. ",
                        .Call("fcnn_ver"), ")\nhttp://fcnn.sourceforge.net/"))
}

#' Fast Compressed Neural Networks for R
#'
#' Provides an interface to kernel routines from the FCNN C++ library.
#' FCNN is based on a completely new Artificial Neural Network
#' representation that offers unmatched efficiency, modularity,
#' and extensibility. FCNN4R provides standard teaching
#' (backpropagation, Rprop, simulated annealing, stochastic gradient)
#' and pruning algorithms (minimum magnitude, Optimal Brain Surgeon),
#' but it is first and foremost an efficient computational engine.
#' Users can easily implement their algorithms by taking advantage
#' of fast gradient computing routines, as well as network
#' reconstruction functionality (removing weights and redundant neurons,
#' reordering inputs, merging networks). Networks can be exported to C functions
#' in order to integrate them into virtually any software solution.
#'
#' @name FCNN4R-package
#'
#' @author Grzegorz Klima <gklima@@users.sourceforge.net>
#'
#' @references
#' G. Klima. \emph{A new approach towards implementing artificial neural networks.}
#' Technical Report, \url{http://fcnn.sourceforge.net/}, 2013.
#'
#' @keywords package
#'
#' @examples
#'
#' # set up the XOR problem inputs and outputs
#' inp <- c(0, 0, 1, 1, 0, 1, 0, 1)
#' dim(inp) <- c(4, 2)
#' outp <- c(0, 1, 1, 0)
#' dim(outp) <- c(4, 1)
#' # create a 2-6-1 network
#' net <- mlp_net(c(2, 6, 1))
#' # set activation function in all layers
#' net <- mlp_set_activation(net, layer = "a", "sigmoid")
#' # randomise weights
#' net <- mlp_rnd_weights(net)
#' # tolerance level
#' tol <- 0.5e-4
#' # teach using Rprop, assign trained network and plot learning history
#' netmse <- mlp_teach_rprop(net, inp, outp, tol_level = tol,
#'                           max_epochs = 500, report_freq = 10)
#' net <- netmse$net
#' plot(netmse$mse, type = 'l')
#' # plot network with weights
#' mlp_plot(net, TRUE)
#' # if the algorithm had converged, prune using Optimal Brain Surgeon and plot
#' if (mlp_mse(net, inp, outp) <= tol) {
#'     net <- mlp_prune_obs(net, inp, outp, tol_level = tol,
#'                          max_reteach_epochs = 500, report = TRUE)[[1]]
#'     mlp_plot(net, TRUE)
#' }
#' # check network output
#' round(mlp_eval(net, inp), digits = 3)
#'
#' @useDynLib FCNN4R
#'
#' @import methods
#' @import Rcpp
#' @importFrom stats runif rnorm
#' @importFrom graphics plot.new plot.window segments text points
#'
NULL
