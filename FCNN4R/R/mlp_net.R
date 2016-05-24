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


#' An S4 class representing Multilayer Perception Network.
#'
#' The \code{mlp_net} class represents the Multilayer Perception Network
#' employing the so-called compressed representation, which was inspired
#' by the Compressed Column Storage familiar from sparse matrix algebra.
#' Although the representation and algorithms working with it are somewhat
#' complicated, the user is provided with a simple and intuitive interface
#' that completely hides the internal workings of the package, which in its
#' large part is written in C++.
#'
#' @aliases mlp_net-class mlp_net-method
#'
#' @slot m_name character string, network name
#' @slot m_layers integer vector, stores the numbers of neurons in layers
#' @slot m_n_pointers integer vector, stores the so-called 'pointers' to neurons
#' @slot m_n_prev integer vector, stores the number of connected neurons in the previous layer
#' @slot m_n_next integer vector, stores the number of connected neurons in the next layer
#' @slot m_w_pointers integer vector, stores the so-called 'pointers' to weights
#' @slot m_w_values numeric vector, values of connection weights and biases
#' @slot m_w_flags logical vector, states (active/inactive) of weights and biases
#' @slot m_w_on integer value, the number of active weights
#' @slot m_af integer vector, activation functions' indices
#' @slot m_af_p numeric vector, activation functions' slope parameters
#'
#' @seealso \code{\link{mlp_net}} for creating objects of this class.
#'
#' @references
#' G. Klima. \emph{A new approach towards implementing artificial neural networks.}
#' Technical Report, \url{http://fcnn.sourceforge.net/}, 2013.
#'
#' @keywords classes
#'
#' @name mlp_net-class
#'
setClass(Class = "mlp_net",
    representation(
        m_name = "character",
        m_layers = "integer",
        m_n_pointers = "integer",
        m_n_prev = "integer",
        m_n_next = "integer",
        m_w_pointers = "integer",
        m_w_values = "numeric",
        m_w_flags = "integer",
        m_w_on = "integer",
        m_af = "integer",
        m_af_p = "numeric"
    ),
    package = "FCNN4R"
)

#' Create objects of \code{mlp_net} class
#'
#' Function used for creating multilayer perceptron networks.
#'
#' @param layers vector providing numbers of neurons in each layer
#' @param name character string, network name (optional)
#'
#' @return Returns an object of \code{mlp_net} class.
#'
#' @seealso \code{\linkS4class{mlp_net}} for details.
#'
#' @examples
#'
#' # create a 2-3-1 network
#' net <- mlp_net(c(2, 3, 1))
#' # randomise weights
#' net <- mlp_rnd_weights(net)
#' # show basic information about the network
#' show(net)
#'
#' @keywords classes
#'
#' @export mlp_net
#'
mlp_net <- function(layers, name = NULL)
{
    layers <- as.integer(layers)
    if (is.null(name)) {
        name <- ""
    } else {
        if (!is.character(name) || (length(name) != 1)) {
            stop("invalid network name")
        }
    }
    cres <- .Call("mlp_construct", layers)
    object <- new("mlp_net",
                  m_name = name,
                  m_layers = layers,
                  m_n_pointers = cres[[1]],
                  m_n_prev = cres[[2]],
                  m_n_next = cres[[3]],
                  m_w_pointers = cres[[4]],
                  m_w_values = cres[[5]],
                  m_w_flags = cres[[6]],
                  m_w_on = cres[[7]],
                  m_af = c(0L, rep(5L, length(layers) - 1)),
                  m_af_p = c(0, rep(.5, length(layers) - 1)))
    return (object)
}


#' Is it an object of \code{mlp_net} class?
#'
#' This function checks whether argument is an object of \code{mlp_net} class.
#'
#' @param x an object to be checked
#'
#' @return Logical value.
#'
#' @keywords classes
#'
#' @export is.mlp_net
#'
is.mlp_net <- function(x)
{
    if (is(x, "mlp_net")) return(TRUE)
    return(FALSE)
}



# #########################################################################
# Displaying networks
# #########################################################################

#' Displaying networks (objects of \code{mlp_net} class)
#'
#' These methods can be used to display objects of \code{mlp_net} class. \code{show}
#' and \code{print} provide short information about network structure and
#' activation functions, \code{summary} gives detailed information about
#' all network connections.
#'
#' @param object an object of \code{mlp_net} class
#' @param x an object of \code{mlp_net} class
#'
#' @name mlp_net-display
#'
#' @aliases show,mlp_net-method
#' @export
#'
setMethod("show", signature(object = "mlp_net"),
function(object)
{
    if (object@m_name != "") {
        cat(paste0("Multilayer perceptron network (", object@m_name, ")\n"))
    } else {
        cat("Multilayer perceptron network\n")
    }
    lays <- object@m_layers
    nlays <- length(lays)
    cat(paste0("Layers: ", lays[1], "(input) - ",
        paste0(lays[2:(nlays - 1)], collapse = " - "),
               " - ", lays[nlays], "(output)\n"))
    cat(paste0("Active weights (connections & biases): ",
               object@m_w_on, " of ", object@m_w_pointers[nlays + 1], "\n"))
    cat("Activation functions:\n")
    for (l in 2:(nlays - 1)) {
        cat(paste0("  layer ", l, " (hidden ", l - 1, "): ",
            mlp_actvfunc2str(object@m_af[l], object@m_af_p[l]), "\n"))
    }
    cat(paste0("  layer ", nlays, " (output): ",
        mlp_actvfunc2str(object@m_af[nlays], object@m_af_p[nlays]), "\n"))
    cat("Weights:\n")
    now <- length(object@m_w_flags)
    truncthresh <- 19
    if (now > truncthresh) {
        trunc <- TRUE
        wg <- as.character(object@m_w_values[1:truncthresh])
        wg[which(object@m_w_flags[1:truncthresh] == 0L)] <- "off"
        wg <- c(wg, "...[truncated]\n")
    } else {
        trunc <- FALSE
        wg <- as.character(object@m_w_values)
        wg[which(object@m_w_flags == 0L)] <- "off"
        wg <- c(wg, "\n")
    }
    cat(wg)
})


#' @rdname mlp_net-display
#' @aliases print,mlp_net-method
#' @export
#'
setMethod("print", signature(x = "mlp_net"),
function(x)
{
    show(x)
})


#' @rdname mlp_net-display
#' @aliases summary,mlp_net-method
#' @export
#'
setMethod("summary", signature(object = "mlp_net"),
function(object)
{
    if (object@m_name != "") {
        cat(paste0("Multilayer perceptron network (", object@m_name, ")\n"))
    } else {
        cat("Multilayer perceptron network\n")
    }
    lays <- object@m_layers
    nlays <- length(lays)
    cat(paste0("Layers: ", lays[1], "(input) - ",
        paste0(lays[2:(nlays - 1)], collapse = " - "),
        " - ", lays[nlays], "(output)\n"))
    cat(paste0("Active weights (connections & biases): ",
               object@m_w_on, " of ", object@m_w_pointers[nlays + 1], "\n"))
    cat("Activation functions:\n")
    for (l in 2:(nlays - 1)) {
        cat(paste0("  layer ", l, " (hidden ", l - 1, "): ",
            mlp_actvfunc2str(object@m_af[l], object@m_af_p[l]), "\n"))
    }
    cat(paste0("  layer ", nlays, " (output): ",
        mlp_actvfunc2str(object@m_af[nlays], object@m_af_p[nlays]), "\n"))
    cat("Weights:\n")
    for (l in 2:(nlays)) {
        if (l < nlays) {
            cat(paste0("  layer ", l, " (hidden layer ", l - 1, "): \n"))
        } else {
            cat(paste0("  layer ", nlays, " (output layer): \n"))
        }
        for (n in (1:lays[l])) {
            cat(paste0("    neuron ", n, ":\n"))
            cat("      bias: ")
            if (!mlp_get_w_st(object, layer = l, nidx = n, nplidx = 0)) {
                cat("off\n")
            } else {
                cat(paste0(mlp_get_w(object, layer = l, nidx = n, nplidx = 0), "\n"))
            }
            for (np in (1:lays[l - 1])) {
                cat(paste0("      conn. to neuron ", np, " in layer ", l - 1, ": "))
                if (!mlp_get_w_st(object, layer = l, nidx = n, nplidx = np)) {
                    cat("off\n")
                } else {
                    cat(paste0(mlp_get_w(object, layer = l, nidx = n, nplidx = np), "\n"))
                }
            }
        }
    }
})




# #########################################################################
# Network names
# #########################################################################


#' Get and set network names
#'
#' The following functions can be used for retrieving and setting network names.
#'
#' @param net an object of \code{mlp_net} class
#' @param name character string with network name
#'
#' @return \code{mlp_get_name} returns character string with network name.
#'
#'         \code{mlp_set_name} returns network (an object of \code{mlp_net}
#'               class) with name set to new value.
#'
#' @name mlp_net-names
#'
#' @export
#'
mlp_get_name <- function(net)
{
    if (!is.mlp_net(net)) {
        stop("expected net argument to be of mlp_net class")
    }
    return(net@m_name)
}



#' @rdname mlp_net-names
#'
#' @export
#'
mlp_set_name <- function(net, name)
{
    if (!is.mlp_net(net)) {
        stop("expected net argument to be of mlp_net class")
    }
    if (!is.character(name) || (length(name) != 1)) {
        stop("invalid network name")
    }
    net@m_name <- name
    return(net)
}



# #########################################################################
# General information about network
# #########################################################################


#' General information about network
#'
#' The following functions return basic information about the network.
#'
#' @param net an object of \code{mlp_net} class
#'
#' @return \code{mlp_get_layers} returns an integer vector with numbers of neurons
#'         in consecutive layers.
#'
#'         \code{mlp_get_no_active_w} returns the number of active weights (connections and biases).
#'
#'         \code{mlp_get_no_w} returns the total number (including inactive) of weights
#'                    (connections and biases).
#'
#' @seealso \code{\link[=mlp_net-class]{mlp_net-class}} for details
#'          on internal network representation.
#'
#' @name mlp_net-general-information
#'
#' @export
#'
mlp_get_layers <- function(net)
{
    if (!is.mlp_net(net)) {
        stop("expected net argument to be of mlp_net class")
    }
    return(net@m_layers)
}

#' @rdname mlp_net-general-information
#'
#' @export
#'
mlp_get_no_active_w <- function(net)
{
    if (!is.mlp_net(net)) {
        stop("expected net argument to be of mlp_net class")
    }
    return(net@m_w_on)
}


#' @rdname mlp_net-general-information
#'
#' @export
#'
mlp_get_no_w <- function(net)
{
    if (!is.mlp_net(net)) {
        stop("expected net argument to be of mlp_net class")
    }
    return(length(net@m_w_flags))
}



# #########################################################################
# Reconstructing network, removing neurons
# #########################################################################

#' Remove redundant neurons in a multilayer perceptron network
#'
#' This function removes redundant neurons from the network, i.e. hidden layers'
#' neurons that are not connected to neurons in the previous layer or the next
#' layer. If a neuron is not connected to neurons in the previous layer but
#' is connected to neurons in the next layer (effectively acts as an additional
#' bias), biases of neurons in the next layer are properly adjusted, therefore,
#' the resulting network behaves just like the initial one.
#'
#' @param net an object of \code{mlp_net} class
#' @param report logical value, if TRUE, information about removed neurons
#'        will be printed on the console (FALSE by default)
#'
#' @return Three-element list. The first element (\code{net}) is the network
#'         (an object of \code{mlp_net} class) with all redundant neurons
#'         removed, the second (\code{ncount}) - the number of neurons removed,
#'         the third (\code{wcount}) - the number of weights removed.
#'
#' @export mlp_rm_neurons
#'
mlp_rm_neurons <- function(net, report = FALSE)
{
    if (!is.mlp_net(net)) {
        stop("expected net argument to be of mlp_net class")
    }
    hst <- net@m_n_pointers[2] + 1
    hen <- net@m_n_pointers[length(net@m_layers)]
    if (all(net@m_n_prev[hst:hen] != 0) && all(net@m_n_next[hst:hen] != 0)) {
        return(list(net = net, ncount = 0, wcount = 0))
    }
    won0 <- net@m_w_on
    rmres <- .Call("mlp_rm_neurons",
                   net@m_layers,
                   net@m_n_pointers,
                   net@m_n_prev,
                   net@m_n_next,
                   net@m_w_pointers,
                   net@m_w_values,
                   net@m_w_flags,
                   net@m_w_on,
                   net@m_af,
                   net@m_af_p,
                   report)
    ret <- new("mlp_net",
               m_name = net@m_name,
               m_layers = rmres[[1]],
               m_n_pointers = rmres[[2]],
               m_n_prev = rmres[[3]],
               m_n_next = rmres[[4]],
               m_w_pointers = rmres[[5]],
               m_w_values = rmres[[6]],
               m_w_flags = rmres[[7]],
               m_w_on = rmres[[8]],
               m_af = net@m_af,
               m_af_p = net@m_af_p)
    return(list(net = ret, ncount = rmres[[9]], wcount = won0 - rmres[[8]]))
}



# #########################################################################
# Manipulating network inputs
# #########################################################################

#' Manipulating network inputs
#'
#' These functions construct new network by removing redundant (i.e. not connected
#' to the next layer) inputs or reordering / expanding network inputs.
#'
#' @param net an object of \code{mlp_net} class
#' @param report logical value, if TRUE, information about removed neurons
#'        will be printed on the console (FALSE by default)
#' @param newnoinputs integer value, determines the number of inputs in the new
#'        network
#' @param inputsmap integer vector, determines the mapping of old inputs into
#'        new ones - the ith value of this vector will be the new index
#'        of ith input
#'
#' @return \code{mlp_rm_input_neurons} returns a two-element list. The first
#'         element (\code{net}) is the network (an object of \code{mlp_net}
#'         class) with all redundant input neurons removed, the second
#'         (\code{ind}) - the indices of input neurons that were not removed.
#'
#'         \code{mlp_expand_reorder_inputs} returns an object of \code{mlp_net}
#'         class.
#'
#' @examples
#'
#' # construct a 2-4-3 network, plot result
#' nn <- mlp_net(c(2, 4, 3))
#' nn <- mlp_rnd_weights(nn)
#' mlp_plot(nn, TRUE)
#' # expand inputs, the new no. of inputs will be 5, with the first input
#' # becoming the 3rd and the second retaining its position, plot result
#' nn <- mlp_expand_reorder_inputs(nn, 5, c(3, 2))
#' mlp_plot(nn, TRUE)
#' # remove redundant neurons (i.e. 1, 4, 5) and plot result
#' nn <- mlp_rm_input_neurons(nn, TRUE)$net
#' mlp_plot(nn, TRUE)
#'
#' @name mlp_net-manipulating-network-inputs
#'
#' @export
#'
mlp_rm_input_neurons <- function(net, report = FALSE)
{
    if (!is.mlp_net(net)) {
        stop("expected net argument to be of mlp_net class")
    }
    if (all(net@m_n_next[1:net@m_layers[1]] != 0)) {
        return(list(net = net, ind = 1:net@m_layers[1]))
    }
    ind <- which(net@m_n_next[1:net@m_layers[1]] != 0)
    rmres <- .Call("mlp_rm_input_neurons",
                   net@m_layers,
                   net@m_n_pointers,
                   net@m_n_prev,
                   net@m_n_next,
                   net@m_w_pointers,
                   net@m_w_values,
                   net@m_w_flags,
                   report)
    ret <- new("mlp_net",
               m_name = net@m_name,
               m_layers = rmres[[1]],
               m_n_pointers = rmres[[2]],
               m_n_prev = rmres[[3]],
               m_n_next = rmres[[4]],
               m_w_pointers = rmres[[5]],
               m_w_values = rmres[[6]],
               m_w_flags = rmres[[7]],
               m_w_on = net@m_w_on,
               m_af = net@m_af,
               m_af_p = net@m_af_p)
    return(list(net = ret, ind = ind))
}


#' @rdname mlp_net-manipulating-network-inputs
#'
#' @export
#'
mlp_expand_reorder_inputs <- function(net, newnoinputs, inputsmap)
{
    if (!is.mlp_net(net)) {
        stop("expected net argument to be of mlp_net class")
    }
    if (!is.numeric(newnoinputs) || length(newnoinputs) != 1) {
        stop("unexpected format of the newnoinputs argument")
    }
    if (!is.numeric(inputsmap) || !is.vector(inputsmap)) {
        stop("unexpected format of the inputsmap argument")
    }
    inputsmap <-  as.integer(inputsmap)
    newnoinputs <- as.integer(newnoinputs)
    rmres <- .Call("mlp_expand_reorder_inputs",
                   net@m_layers,
                   net@m_n_pointers,
                   net@m_n_prev,
                   net@m_n_next,
                   net@m_w_pointers,
                   net@m_w_values,
                   net@m_w_flags,
                   newnoinputs,
                   inputsmap)
    ret <- new("mlp_net",
               m_name = net@m_name,
               m_layers = rmres[[1]],
               m_n_pointers = rmres[[2]],
               m_n_prev = rmres[[3]],
               m_n_next = rmres[[4]],
               m_w_pointers = rmres[[5]],
               m_w_values = rmres[[6]],
               m_w_flags = rmres[[7]],
               m_w_on = net@m_w_on,
               m_af = net@m_af,
               m_af_p = net@m_af_p)
    return(ret)
}







# #########################################################################
# Combining two networks into one
# #########################################################################

#' Combining two networks into one
#'
#' These functions construct new network by merging two networks
#' (they must have the same number of layers) or by connecting
#' one network outputs to another network inputs (the numbers of output
#' and input neurons must agree). These functions may be used in constructing
#' deep learning networks or constructing networks with some special topologies.
#'
#' @param net1 an object of \code{mlp_net} class
#' @param net2 an object of \code{mlp_net} class
#' @param same_inputs logical, if TRUE both merged networks are assumed to take
#'          the same inputs (they share the input layer), default is FALSE
#'
#' @return Both functions return an object of \code{mlp_net} class.
#'
#' @examples
#'
#' # create two 2-2-2 networks with random weights and plot them
#' net1 <- mlp_net(c(2, 2, 2))
#' net1 <- mlp_rnd_weights(net1)
#' mlp_plot(net1, TRUE)
#' net2 <- mlp_net(c(2, 2, 2))
#' net2 <- mlp_rnd_weights(net2)
#' mlp_plot(net2, TRUE)
#' # create a 4-3-2 network with random weights and plot it
#' net3 <- mlp_net(c(4, 3, 2))
#' net3 <- mlp_rnd_weights(net3)
#' mlp_plot(net3, TRUE)
#' # construct new network using net1, net2, and net3 and plot it
#' net4 <- mlp_stack(mlp_merge(net1, net2), net3)
#' mlp_plot(net4, TRUE)
#'
#' @name mlp_net-combining-two-networks
#'
#' @export
#'
mlp_merge <- function(net1, net2, same_inputs = FALSE)
{
    if (!is.mlp_net(net1)) {
        stop("expected net1 argument to be of mlp_net class")
    }
    if (!is.mlp_net(net2)) {
        stop("expected net2 argument to be of mlp_net class")
    }
    if (!is.logical(same_inputs)) {
        stop("expected logical argument")
    }
    if (length(net1@m_layers) != length(net2@m_layers)) {
        stop("different numbers of layers in networks");
    }
    if (any(net1@m_af != net2@m_af) || any(net1@m_af_p != net2@m_af_p)) {
        stop("activation functions in networks disagree");
    }
    res <- .Call("mlp_merge",
                 net1@m_layers,
                 net1@m_w_pointers,
                 net1@m_w_values,
                 net1@m_w_flags,
                 net2@m_layers,
                 net2@m_w_pointers,
                 net2@m_w_values,
                 net2@m_w_flags,
                 same_inputs)
    net <- new("mlp_net",
               m_name = "",
               m_layers = res[[1]],
               m_n_pointers = res[[2]],
               m_n_prev = res[[3]],
               m_n_next = res[[4]],
               m_w_pointers = res[[5]],
               m_w_values = res[[6]],
               m_w_flags = res[[7]],
               m_w_on = res[[8]],
               m_af = net1@m_af,
               m_af_p = net1@m_af_p)
    return(net)
}


#' @rdname mlp_net-combining-two-networks
#'
#' @export
#'
mlp_stack <- function(net1, net2)
{
    if (!is.mlp_net(net1)) {
        stop("expected net1 argument to be of mlp_net class")
    }
    if (!is.mlp_net(net2)) {
        stop("expected net2 argument to be of mlp_net class")
    }
    res <- .Call("mlp_stack",
                 net1@m_layers,
                 net1@m_w_pointers,
                 net1@m_w_values,
                 net1@m_w_flags,
                 net2@m_layers,
                 net2@m_w_pointers,
                 net2@m_w_values,
                 net2@m_w_flags)
    net <- new("mlp_net",
               m_name = "",
               m_layers = res[[1]],
               m_n_pointers = res[[2]],
               m_n_prev = res[[3]],
               m_n_next = res[[4]],
               m_w_pointers = res[[5]],
               m_w_values = res[[6]],
               m_w_flags = res[[7]],
               m_w_on = res[[8]],
               m_af = c(net1@m_af, net2@m_af[2:length(net2@m_af)]),
               m_af_p = c(net1@m_af_p, net2@m_af_p[2:length(net2@m_af_p)]))
    return(net)
}



# #########################################################################
# Importing and exporting networks
# #########################################################################

#' Export and import multilayer perceptron network to/from a text file
#' in FCNN format
#'
#' These functions can be used to export and import multilayer perceptron
#' network to/from a text file in FCNN format.
#'
#' Files are organised as follows:
#' \itemize{
#'  \item the first comment (beginning with \code{#}) is treated as network
#'  information (name) string,
#'  \item all other comments are ignored,
#'  \item network structure is represented by five block of numbers:
#'     \itemize{
#'      \item the first line determines numbers of neurons in consecutive layers,
#'      \item the second block of 0's and 1's determines which weights are turned off/on,
#'      \item the third block contains active weights' values,
#'      \item the last block determines hidden and output layers' activation functions
#'            and their slope parameters - each line contains 2 numbers: the function index
#'            and its slope parameter.
#'      }
#'  }
#'
#' @param fname character string with the filename
#' @param net an object of \code{mlp_net} class
#'
#' @return \code{mlp_export_fcnn} returns logical value, TRUE if export was successful,
#'         FALSE otherwise.
#'
#'         \code{mlp_import_fcnn} returns an object of \code{mlp_net} class or NULL,
#'         if import failed.
#'
#' @seealso \code{\linkS4class{mlp_net}} for network representation details.
#'
#' @examples
#'
#' # create a 2-3-1 network
#' net <- mlp_net(c(2, 3, 1))
#' # randomise weights
#' net <- mlp_rnd_weights(net)
#' # show the network
#' show(net)
#' # export network
#' mlp_export_fcnn("test.net", net)
#' # show the output file
#' file.show("test.net")
#' # import network
#' net2 <- mlp_import_fcnn("test.net")
#' # show the imported network
#' show(net2)
#'
#' @name mlp_net-export-import
#'
#' @export
#'
mlp_export_fcnn <- function(fname, net)
{
    if (!is.character(fname) || (length(fname) != 1)) {
        stop("invalid filename")
    }
    if (!is.mlp_net(net)) {
        stop("expected net argument to be of mlp_net class")
    }
    return(.Call("mlp_export", fname,
                 net@m_name, net@m_layers, net@m_w_values, net@m_w_flags,
                 net@m_af, net@m_af_p))
}

#' @rdname mlp_net-export-import
#'
#' @export
#'
mlp_import_fcnn <- function(fname)
{
    if (!is.character(fname) || (length(fname) != 1)) {
        stop("invalid filename")
    }
    impres <- .Call("mlp_import", fname)
    if (is.null(impres)) {
        return(impres)
    }
    object <- new("mlp_net",
                  m_name = impres[[1]],
                  m_layers = impres[[2]],
                  m_n_pointers = impres[[3]],
                  m_n_prev = impres[[4]],
                  m_n_next = impres[[5]],
                  m_w_pointers = impres[[6]],
                  m_w_values = impres[[7]],
                  m_w_flags = impres[[8]],
                  m_w_on = impres[[9]],
                  m_af = impres[[10]],
                  m_af_p = impres[[11]])
    return(object)
}



# #########################################################################
# Exporting networks to C
# #########################################################################

#' Export multilayer perceptron network to a C function
#'
#' This function exports multilayer perceptron network to a C function
#' with optional affine input and output transformations: Ax+b for inputs
#' and Cx+d for outputs.
#'
#' @param fname character string with the filename
#' @param net an object of \code{mlp_net} class
#' @param with_bp logical, should backpropagation code for online learning
#'                be exported?
#' @param A numeric matrix (optional), input linear transformation
#' @param b numeric vector (optional), input translation
#' @param C numeric matrix (optional), output linear transformation
#' @param d numeric vector (optional), output translation
#'
#' @return Logical value, TRUE if export was successful, FALSE otherwise.
#'
#' @examples
#'
#' # create a 2-3-1 network
#' net <- mlp_net(c(2, 3, 1))
#' # randomise weights
#' net <- mlp_rnd_weights(net)
#' # show the network
#' show(net)
#' # export network to a C function
#' mlp_export_C("test.c", net)
#' # show the output file
#' file.show("test.c")
#'
#' @export mlp_export_C
#'
mlp_export_C <- function(fname, net, with_bp = FALSE, A = NULL, b = NULL, C = NULL, d = NULL)
{
    if (!is.character(fname) || (length(fname) != 1)) {
        stop("invalid filename")
    }
    if (!is.mlp_net(net)) {
        stop("expected net argument to be of mlp_net class")
    }
    if (!is.logical(with_bp) || (length(with_bp) != 1)) {
        stop("expected with_bp argument to be logical value")
    }

    if (!is.null(A)) {
        if (!is.numeric(A) || (length(dim(A)) != 2)) {
            stop("A must be a numeric matrix")
        }
        if (!all(dim(A) == net@m_layers[1]))
            stop("invalid sizes of matrix A")
    }
    if (!is.null(b)) {
        if (!is.numeric(b)) {
            stop("b must be a numeric vector")
        }
        if (length(b) != net@m_layers[1])
            stop("invalid size of vector b")
    }
    if ((!is.null(A) || !is.null(b)) && ((is.null(A) || is.null(b)))) {
        stop("incomplete input transformation provided")
    }
    if (!is.null(C)) {
        if (!is.numeric(C) || (length(dim(C)) != 2)) {
            stop("C must be a numeric matrix")
        }
        if (!all(dim(C) == net@m_layers[length(net@m_layers)]))
            stop("invalid sizes of matrix C")
    }
    if (!is.null(d)) {
        if (!is.numeric(d)) {
            stop("d must be a numeric vector")
        }
        if (length(d) != net@m_layers[length(net@m_layers)])
            stop("invalid size of vector d")
    }
    if ((!is.null(C) || !is.null(d)) && ((is.null(C) || is.null(d)))) {
        stop("incomplete output transformation provided")
    }
    if (with_bp && !is.null(C)) {
        tryCatch({E <- solve(C);}, error = function(e) stop("output transformation is not invertible"))
        f <- -E %*% d
    } else {
        E <- NULL
        f <- NULL
    }

    return(.Call("mlp_export_C", fname,
                 net@m_name, net@m_layers, net@m_n_pointers,
                 net@m_w_values, net@m_w_flags, net@m_w_on,
                 net@m_af, net@m_af_p, with_bp,
                 A, b, C, d, E, f))
}



# #########################################################################
# Activation functions
# #########################################################################


#' Return character string representing activation function
#'
#' @param idx activation function index
#' @param slope activation function index slope parameter
#'
#' @return This function returns character string representing activation function.
#'
#' @keywords internal
#'
mlp_actvfunc2str <- function(idx, slope)
{
    strng <- .Call("actvfuncstr", idx)
    if (idx > 2) strng <- paste0(strng, " with s = ", slope)
    return (strng)
}



#' Set network activation functions
#'
#' This function sets activation function (and its slope parameter)
#' for neurons in the hidden layers and in the output layer.
#'
#' @param net an object of \code{mlp_net} class
#' @param layer integer vector or character value, index (indices) of layer(s)
#'              whose activation function will be changed or character:
#'              "a" denotes all layers, "h" - hidden layer(s), "o" - the output layer
#' @param activation character string, activation function name, admissible
#'                   options are: "threshold", "sym_threshold", "linear",
#'                   "sigmoid", "sym_sigmoid" (and "tanh"), "sigmoid_approx",
#'                   and "sym_sigmoid_approx"
#' @param slope numeric value, activation function slope parameter, if 0
#'              the default parameter value is chosen for each activation function
#'
#' @return This function returns network (an object of \code{mlp_net} class)
#'         with activation function set.
#'
#' @export mlp_set_activation
#'
mlp_set_activation <- function(net,
                               layer,
                               activation = c("threshold", "sym_threshold", "linear",
                                              "sigmoid", "sym_sigmoid", "tanh",
                                              "sigmoid_approx", "sym_sigmoid_approx"),
                               slope = 0)
{
    if (!is.mlp_net(net)) {
        stop("expected net argument to be of mlp_net class")
    }
    if (is.character(layer)) {
        if (layer == "a") {
            layer <- 2:length(net@m_layers)
        } else if (layer == "h") {
            layer <- 2:(length(net@m_layers) - 1)
        } else if (layer == "o") {
            layer <- length(net@m_layers)
        } else {
            stop("invalid layer argument - character should be \"a\", \"h\", or \"o\"")
        }
    } else if (is.numeric(layer)) {
        if (!all(layer %in% 2:length(net@m_layers))) {
            stop("invalid layer index")
        }
    } else {
        stop("expected layer argument to be integer or character value")
    }

    if ((length(activation) != 1) && (length(activation) != length(layer))) {
        stop("incompatible lengths of layer and activation arguments")
    }
    if ((length(activation) == 1) && (length(layer) != 1)) {
        activation <- rep(activation, length(layer))
    }

    afi <- rep(0L, length(activation))
    afi[which(activation == "threshold")] <- 1L
    afi[which(activation == "sym_threshold")] <- 2L
    afi[which(activation == "linear")] <- 3L
    afi[which(activation == "sigmoid")] <- 4L
    afi[which(activation == "sym_sigmoid")] <- 5L
    afi[which(activation == "tanh")] <- 5L
    afi[which(activation == "sigmoid_approx")] <- 6L
    afi[which(activation == "sym_sigmoid_approx")] <- 7L

    if (any(afi == 0L)) {
        stop("invalid activation function name")
    }

    if (!is.numeric(slope) || any(!is.finite(slope)) || any(slope < 0)) {
        stop("invalid slope parameter")
    }
    if ((length(slope) != 1) && (length(slope) != length(layer))) {
        stop("incompatible lengths of layer and slope arguments")
    }
    if ((length(slope) == 1) && (length(layer) != 1)) {
        slope <- rep(slope, length(layer))
    }
    slope[which((slope == 0) & (afi %in% 1:3))] <- 1
    slope[which((slope == 0) & (afi %in% 4:7))] <- 0.5
    net@m_af[layer] <- afi
    net@m_af_p[layer] <- slope
    return(net)
}



# #########################################################################
# Weights indexing
# #########################################################################


#' Check validity of weight index
#'
#' @param net an object of \code{mlp_net} class
#' @param idx integer value (vector), weight absolute index
#' @param layer integer value (vector), layer index
#' @param nidx integer value (vector), neuron index
#' @param nplidx integer value (vector), index of the neuron in the previous
#'               layer determining connection from neuron \code{nidx}
#'               in \code{layer}, 0 denotes bias of neuron \code{nidx} in \code{layer}
#'
#' @return This function does not return.
#'
#' @keywords internal
#'
mlp_check_w <- function(net, idx = NULL, layer = NULL, nidx = NULL, nplidx = NULL)
{
    if (!is.mlp_net(net)) {
        stop("expected net argument to be of mlp_net class")
    }
    if (!is.null(idx)) {
        if (any(idx < 1) || any(idx > length(net@m_w_values))) {
            stop("invalid weight index")
        }
        if (!is.null(layer) || !is.null(nidx) || !is.null(nplidx)) {
            stop("weight idx already provided, other arguments should not be set")
        }
    } else {
        if (is.null(layer) || is.null(nidx) || is.null(nplidx)) {
            stop("weight idx not given, 3 arguments (layer, nidx, nplidx) required")
        }
        if (any(layer < 2) || any(layer > length(net@m_layers))) {
            stop("invalid layer")
        }
        if (any(nidx < 1) || any(nidx > net@m_layers[layer])) {
            stop("invalid neuron index (nidx)")
        }
        if (any(nplidx < 0) || any(nplidx > net@m_layers[layer - 1])) {
            stop("invalid previous layer neuron index (nplidx)")
        }
    }
}


#' Retrieving absolute weight index
#'
#' In some situations absolute weight index (i.e. index within all weights
#' including inactive ones) needs to be computed based on information
#' about connected neurons' indices or weight index within actives ones.
#' The latter functionality is especially useful in implementation of pruning
#' algorithms.
#'
#' @param net an object of \code{mlp_net} class
#' @param layer integer value (vector), layer index
#' @param nidx integer value (vector), neuron index
#' @param nplidx integer value (vector), index of the neuron in the previous
#'               layer determining connection from neuron \code{nidx}
#'               in \code{layer}, 0 denotes bias of neuron \code{nidx} in \code{layer}
#' @param idx integer value (vector), weight index (indices) within active ones
#'
#' @return Absolute weight index.
#'
#' @name mlp_net-absolute-weight-indices
#'
#' @export
#'
mlp_get_w_idx <- function(net, layer, nidx, nplidx)
{
    mlp_check_w(net, idx = NULL, layer = layer, nidx = nidx, nplidx = nplidx)
    idx <- net@m_w_pointers[layer] +
           (nidx - 1) * (net@m_layers[layer - 1] + 1) + nplidx + 1
    return(idx)
}


#' @name mlp_net-absolute-weight-indices
#'
#' @export
#'
mlp_get_w_abs_idx <- function(net, idx)
{
    if (!is.mlp_net(net)) {
        stop("expected net argument to be of mlp_net class")
    }
    if (any(idx < 1) || any(idx > net@m_w_on)) {
        stop("invalid active weight index")
    }
    return(.Call("mlp_get_abs_w_idx", net@m_w_flags, as.integer(idx)))
}



# #########################################################################
# Individual weight access
# #########################################################################

#' Setting and retrieving status (on/off) and value of individual weight(s)
#'
#' The following functions can be used to access individual weight(s), i.e. set or
#' retrieve status(es) (on/off) and value(s).
#'
#' @param net an object of \code{mlp_net} class
#' @param on logical value (vector), should the weight be set on or off?
#' @param val numeric value (vector), connection (or bias) value to be set
#' @param idx integer value (vector), weight absolute index
#' @param layer integer value (vector), layer index
#' @param nidx integer value (vector), neuron index
#' @param nplidx integer value (vector), index of the neuron in the previous
#'               layer determining connection from neuron \code{nidx}
#'               in \code{layer}, 0 denotes bias of neuron \code{nidx} in \code{layer}
#'
#' @return \code{mlp_set_w_st} returns network (an object of \code{mlp_net} class)
#'                    with state(s) (on/off) of selected weight(s) set.
#'
#'         \code{mlp_set_w} returns network (an object of \code{mlp_net} class)
#'                    with value(s) of selected weight(s) set.
#'
#'         \code{mlp_get_w_st} returns logical value (vector), TRUE if connection/bias is active,
#'                    FALSE otherwise.
#'
#'         \code{mlp_get_w} returns numeric value (vector), selected weight value(s).
#'
#' @name mlp_net-accessing-individual-weights
#'
#' @export
#'
mlp_set_w_st <- function(net, on, idx = NULL, layer = NULL, nidx = NULL, nplidx = NULL)
{
    mlp_check_w(net, idx, layer, nidx, nplidx)
    if (is.null(idx)) {
        idx <- mlp_get_w_idx(net, layer, nidx, nplidx)
    }
    N <- length(idx)
    if (length(on) == 1) {
        on <- rep(on, N)
    } else if (length(on) != N) {
        stop("nonconformant lengths of inputs")
    }
    output <- .C("mlp_set_active",
                 net@m_layers, net@m_n_pointers,
                 net@m_n_prev, net@m_n_next,
                 net@m_w_pointers, net@m_w_values,
                 net@m_w_flags, net@m_w_on,
                 as.integer(idx), as.integer(on), as.integer(N))
    net@m_n_prev <- output[[3]]
    net@m_n_next <- output[[4]]
    net@m_w_values <- output[[6]]
    net@m_w_flags <- output[[7]]
    net@m_w_on <- output[[8]]
    return(net)
}


#' @rdname mlp_net-accessing-individual-weights
#'
#' @export
#'
mlp_set_w <- function(net, val, idx = NULL, layer = NULL, nidx = NULL, nplidx = NULL)
{
    mlp_check_w(net, idx, layer, nidx, nplidx)
    if (is.null(idx)) {
        idx <- mlp_get_w_idx(net, layer, nidx, nplidx)
    }
    if (any(net@m_w_flags[idx] == 0L)) {
        if (length(idx) == 1) {
            stop("selected weight is off")
        } else {
            stop("at least one selected weight is off")
        }
    }
    net@m_w_values[idx] <- val
    return(net)
}


#' @rdname mlp_net-accessing-individual-weights
#'
#' @export
#'
mlp_get_w_st <- function(net, idx = NULL, layer = NULL, nidx = NULL, nplidx = NULL)
{
    mlp_check_w(net, idx, layer, nidx, nplidx)
    if (is.null(idx)) {
        idx <- mlp_get_w_idx(net, layer, nidx, nplidx)
    }
    return(net@m_w_flags[idx] != 0L)
}


#' @rdname mlp_net-accessing-individual-weights
#'
#' @export
#'
mlp_get_w <- function(net, idx = NULL, layer = NULL, nidx = NULL, nplidx = NULL)
{
    mlp_check_w(net, idx, layer, nidx, nplidx)
    if (is.null(idx)) {
        idx <- mlp_get_w_idx(net, layer, nidx, nplidx)
    }
    if (any(net@m_w_flags[idx] == 0L)) {
        if (length(idx) == 1) {
            stop("selected weight is off")
        } else {
            stop("at least one selected weight is off")
        }
    }
    return(net@m_w_values[idx])
}


# #########################################################################
# Working with weights' vector
# #########################################################################

#' This function sets network weights to random values drawn from uniform
#' distribution.
#'
#' @param net an object of \code{mlp_net} class
#' @param a numeric value, values will be drawn from uniform distribution
#'        on [-a, a] (by default a = 0.2)
#'
#' @return Network (an object of \code{mlp_net} class) with randomised weights.
#'
#' @export mlp_rnd_weights
#'
mlp_rnd_weights <- function(net, a = 0.2)
{
    if (!is.mlp_net(net)) {
        stop("expected net argument to be of mlp_net class")
    }
    ind <- which(net@m_w_flags != 0L)
    weights <- runif(length(ind), min = -a, max = a)
    net@m_w_values[ind] <- weights
    return(net)
}


#' Set and retrieve (active) weights' values
#'
#' One of FCNN's design objectives (and main advantages) is the complete separation
#' of teaching (and pruning) algorithms from internal network structure workings.
#' This goal is achieved through fast access to (active) weights vector facilitated
#' by FCNN's `compressed' network representation. The following two functions
#' allow users to efficiently retrieve and set network (active) weights vector.
#'
#' @param net an object of \code{mlp_net} class
#' @param weights numeric vector of new active weights' values
#'
#' @return \code{mlp_set_weights} returns network (an object of \code{mlp_net}
#'         class) with active weights set to given values.
#'
#'         \code{mlp_set_weights} returns numeric vector of active weights' values.
#'
#' @name mlp_net-weights-access
#'
#' @export
#'
mlp_set_weights <- function(net, weights)
{
    if (!is.mlp_net(net)) {
        stop("expected net argument to be of mlp_net class")
    }
    if (length(weights) != net@m_w_on) {
        stop("invalid size of active weights vector")
    }
    ind <- which(net@m_w_flags != 0L)
    net@m_w_values[ind] <- weights
    return(net)
}


#' @rdname mlp_net-weights-access
#'
#' @export
#'
mlp_get_weights <- function(net)
{
    if (!is.mlp_net(net)) {
        stop("expected net argument to be of mlp_net class")
    }
    ind <- which(net@m_w_flags != 0L)
    return(net@m_w_values[ind])
}





# #########################################################################
# Evaluation, MSE, and gradients
# #########################################################################


#' Check validity of inputs and outputs
#'
#' @param net an object of \code{mlp_net} class
#' @param input numeric matrix, each row corresponds to one input vector,
#'        the number of columns must be equal to the number of neurons
#'        in the network input layer
#' @param output numeric matrix with rows corresponding to expected outputs,
#'        the number of columns must be equal to the number of neurons
#'        in the network output layer, the number of rows must be equal to the number
#'        of input rows
#' @param i data row index
#'
#' @return This function does not return.
#'
#' @keywords internal
#'
mlp_check_inout <- function(net, input, output = NULL, i = NULL)
{
    if (!is.mlp_net(net)) {
        stop("expected net argument to be of mlp_net class")
    }
    if (!is.numeric(input)) {
        stop("invalid input, expected numeric matrix")
    }
    di <- dim(input)
    if (length(di) != 2) {
        stop("invalid input, expected numeric matrix")
    }
    if (di[1] == 0) {
        stop("input data must have at least one row")
    }
    if (di[2] != net@m_layers[1]) {
        stop("number of input columns does not match the number of input neurons")
    }
    if (!is.null(i)) {
        if ((length(i) != 1) || (!is.integer(i) && !is.numeric(i))) {
            stop("invalid type of row index")
        }
        if ((i < 1) || (i > di[1])) {
            stop("invalid row index")
        }
    }
    if (!is.null(output)) {
        if (!is.numeric(output)) {
            stop("invalid output, expected numeric matrix")
        }
        do <- dim(output)
        if (length(do) != 2) {
            stop("invalid output, expected numeric matrix")
        }
        if (do[1] != di[1]) {
            stop("no. of output rows and no. of input rows disagree")
        }
        if (do[2] != net@m_layers[length(net@m_layers)]) {
            stop("number of output columns does not match the number of output neurons")
        }
    }
}



#' Evaluation
#'
#' Evaluate network output.
#'
#' @param net an object of \code{mlp_net} class
#' @param input numeric matrix, each row corresponds to one input vector,
#'        the number of columns must be equal to the number of neurons
#'        in the network input layer
#'
#' @return Numeric matrix with rows representing network outputs corresponding
#'         to input rows.
#'
#' @export mlp_eval
#'
mlp_eval <- function(net, input)
{
    mlp_check_inout(net, input)
    lays <- net@m_layers
    nout <- lays[length(lays)]
    nrows <- dim(input)[1]
    output <- matrix(0, nrow = nrows, ncol = nout)
    output <- .C("mlp_eval",
                 as.integer(lays), as.integer(length(lays)), as.integer(net@m_n_pointers),
                 as.numeric(net@m_w_values),
                 as.integer(net@m_af), as.numeric(net@m_af_p),
                 as.integer(nrows), as.numeric(input), res = as.numeric(output))$res
    dim(output) <- c(nrows, nout)
    return(output)
}



#' Computing mean squared error, its gradient, and output derivatives
#'
#' The functions use fast FCNN kernel routines and are intended for implementing
#' teaching and pruning algorithms.
#'
#' \code{mlp_mse} returns the mean squared error (MSE). MSE is understood
#' as half of the squared error averaged over all outputs and data records.
#'
#' \code{mlp_grad} computes the gradient of MSE w.r.t. network weights.
#' This function is useful when implementing batch teaching algorithms.
#'
#' \code{mlp_gradi} computes the gradient of MSE w.r.t. network weights at the \code{i}th
#' data record. This is normalised by the number of outputs only,
#' the average over all rows (all i) returns the same as \code{grad(input, output)}.
#' This function is useful for implementing on-line teaching algorithms.
#'
#' \code{mlp_gradij} computes gradients of network outputs,
#' i.e the derivatives of outputs w.r.t. active weights, at given data row.
#' The derivatives of outputs are placed in subsequent columns of the returned
#' matrix. Scaled by the output errors and averaged they give the same
#' as \code{gradi(input, output, i)}. This function is useful in implementing
#' teaching algorithms using second order corrections and Optimal Brain Surgeon
#' pruning algorithm.
#'
#' \code{mlp_jacob} computes the Jacobian of network outputs, i.e the derivatives
#' of outputs w.r.t. inputs, at given data row.
#' The derivatives of outputs are placed in subsequent columns of the returned
#' matrix.
#'
#' @param net an object of \code{mlp_net} class
#' @param input numeric matrix, each row corresponds to one input vector,
#'        the number of columns must be equal to the number of neurons
#'        in the network input layer
#' @param output numeric matrix with rows corresponding to expected outputs,
#'        the number of columns must be equal to the number of neurons
#'        in the network output layer, the number of rows must be equal to the number
#'        of input rows
#' @param i data row index
#'
#' @return \code{mlp_mse} returns mean squared error (numeric value).
#'
#' \code{mlp_grad} returns two-element lists with the first
#' field (\code{grad}) containing numeric vector with gradient and the second
#' (\code{mse}) - the mean squared error.
#'
#' \code{mlp_gradi} returns numeric vector with gradient.
#'
#' \code{mlp_gradij} returns numeric matrix with gradients of outputs in
#' consecutive columns.
#'
#' \code{mlp_jacob} returns numeric matrix with derivatives of outputs in
#' consecutive columns.
#'
#' @name mlp_net-MSE-gradients
#'
#' @export mlp_mse
#'
mlp_mse <- function(net, input, output)
{
    mlp_check_inout(net, input, output = output)
    lays <- net@m_layers
    nout <- lays[length(lays)]
    nrows <- dim(input)[1]
    mse <- 0
    mse <- .C("mlp_mse",
              as.integer(lays), as.integer(length(lays)), as.integer(net@m_n_pointers),
              as.numeric(net@m_w_values),
              as.integer(net@m_af), as.numeric(net@m_af_p),
              as.integer(nrows), as.numeric(input),
              as.numeric(output), res = mse)$res
    return(mse)
}



#' @rdname mlp_net-MSE-gradients
#'
#' @export
#'
mlp_grad <- function(net, input, output)
{
    mlp_check_inout(net, input, output = output)
    lays <- net@m_layers
    nout <- lays[length(lays)]
    nrows <- dim(input)[1]
    grad <- numeric(length = net@m_w_on + 1)
    grad <- .C("mlp_grad",
               lays, length(lays), net@m_n_pointers,
               net@m_w_pointers, net@m_w_flags, net@m_w_values,
               as.integer(net@m_af), as.numeric(net@m_af_p),
               as.integer(nrows), as.numeric(input), as.numeric(output),
               res = grad)$res
    mse <- grad[1]
    grad <- grad[2:length(grad)]
    return(list(grad = grad, mse = mse))
}



#' @rdname mlp_net-MSE-gradients
#'
#' @export
#'
mlp_gradi <- function(net, input, output, i)
{
    mlp_check_inout(net, input, output = output, i = i)
    lays <- net@m_layers
    nout <- lays[length(lays)]
    nrows <- dim(input)[1]
    grad <- numeric(length = net@m_w_on)
    grad <- .C("mlp_gradi",
               lays, length(lays), net@m_n_pointers,
               net@m_w_pointers, net@m_w_flags, net@m_w_values,
               as.integer(net@m_af), as.numeric(net@m_af_p),
               as.integer(nrows), as.integer(i), as.numeric(input), as.numeric(output),
               res = grad)$res
    return(grad)
}



#' @rdname mlp_net-MSE-gradients
#'
#' @export
#'
mlp_gradij <- function(net, input, i)
{
    mlp_check_inout(net, input, i = i)
    lays <- net@m_layers
    nout <- lays[length(lays)]
    nrows <- dim(input)[1]
    grad <- matrix(0, nrow = net@m_w_on, ncol = nout)
    grad <- .C("mlp_gradij",
               lays, length(lays), net@m_n_pointers,
               net@m_w_pointers, net@m_w_flags,
               net@m_w_values, net@m_w_on,
               as.integer(net@m_af), as.numeric(net@m_af_p),
               as.integer(nrows), as.integer(i), as.numeric(input),
               res = grad)$res
    dim(grad) <- c(net@m_w_on, nout)
    return(grad)
}



#' @rdname mlp_net-MSE-gradients
#'
#' @export
#'
mlp_jacob <- function(net, input, i)
{
    mlp_check_inout(net, input, i = i)
    lays <- net@m_layers
    nin <- lays[1]
    nout <- lays[length(lays)]
    nrows <- dim(input)[1]
    jacob <- matrix(0, nrow = nin, ncol = nout)
    jacob <- .C("mlp_jacob",
                lays, length(lays), net@m_n_pointers,
                net@m_w_pointers, net@m_w_flags,
                net@m_w_values, net@m_w_on,
                as.integer(net@m_af), as.numeric(net@m_af_p),
                as.integer(nrows), as.integer(i), as.numeric(input),
                res = jacob)$res
    dim(jacob) <- c(nin, nout)
    return(jacob)
}






