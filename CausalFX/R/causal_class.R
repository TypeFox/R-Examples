############################################################################################
# causal_class.R
#
# Basic definitions for causal effect learning algorithms.
#
# Code by
#
#  - Ricardo Silva (ricardo@stats.ucl.ac.uk)
#  - Robin Evans (robin.evans@stats.ox.ac.uk)
#
# Current version: 28/04/2015
# First version: 31/03/2014

#' @title Creates a CausalFX Problem Instance
#' 
#' @description 
#' Set up an object describing a causal inference problem of finding the
#' average causal effect of some treatment on some outcome. Currently, only
#' binary data is supported. The problem specification also allows the specification
#' of a synthetic model, for simulation studies.
#'
#'@param x the index of the treatment variable.
#'@param y the index of the outcome variable.
#'@param latent_idx an array with the indices of variables which should be considered latent
#'@param dat a matrix of binary data, can be ignored if a model is provided.
#'@param g a binary matrix encoding a causal graph, where g[i, j] == 1 if 
#'       a directed edge from vertex \eqn{j} to \eqn{i} should exist, 0 otherwise. This is only required
#'       if a ground truth model exists.
#'@param model if \code{g} is specified, this needs to be specified too. This argument should be 
#'             a list of conditional probability tables, each encoding the conditional probability of 
#'             each vertex in \code{g} given its parents. Entry \code{model[[i]]} is 
#'             an array of non-negative numbers, describing the probability of random variable/vertex
#'             \eqn{i} being equal to 1. In particular, \code{model[[i]][j]} is the conditional
#'             probability of this event given that the parents of \eqn{i} are in state \eqn{j}.
#'             States are indexed as follows. If \eqn{S} is the binary string corresponding to the 
#'             binary values of the parents of \eqn{i} in \code{g}, sorted by their index, then
#'             \eqn{j} is given by \eqn{1 + bin2dec(S)}, where \eqn{bin2dec} is the transformation of
#'             a binary string into a decimal number.
#'@param num_v_max the maximum dimensionality in which the joint distribution implied by a model is pre-computed.
#'                 Having this pre-computed can speed up some computations for methods that use the provided
#'                 ground truth model. Because the space required to store a joint distribution
#'                 grows exponentially with the dimensionality, this quantity cannot be too large.
#'
#' @return A \code{cfx} object, which contains the following fields:
#'  \item{\code{X_idx}}{the index of the treatment variable in the data/graph.}
#'  \item{\code{Y_idx}}{the index of the outcome variable in the data/graph.}
#'  \item{\code{latent_idx}}{the array of latent variable indices given as input.}
#'  \item{\code{data}}{the data given as input.}
#'  \item{\code{graph}}{the graph given as input.}
#'  \item{\code{varnames}}{an array of strings with the names of the variables, as given by \code{data}. 
#'                         If \code{data} has no column names or it is not provided, this is given a default value, where
#'                         variable \eqn{i} is assigned the name "\code{X}".}
#'  \item{\code{model}}{the model given as input.}
#'  \item{\code{ancestrals}}{a list of arrays (if \code{g} is provided), where \code{ancestrals[[i]]} is the array of the indices of 
#'                           the ancestrals of \eqn{i} in \code{g}, excluding \eqn{i} itself.}
#'  \item{\code{probs}}{a multidimensional array (if \code{g} is provided) of the same dimensionality as
#'                      \code{g}, where each entry corresponds to the probability of that particular assignment of variable values.
#'                      This is \code{NULL} if the dimensionality of \code{g} is greater than \code{num_v_max}.}
#'
#'@export

cfx <- function(x, y, latent_idx = NULL, dat = NULL, g = NULL, model = NULL, num_v_max = 20) {

  if (is.null(latent_idx)) latent_idx <- c()
  
  out <- list(X_idx = x, Y_idx = y, data = dat, graph = g, model = model, latent_idx = latent_idx)
  
  if (!is.null(dat)) {
    num_v <- ncol(dat)
  } else if (!is.null(g)) {
    num_v <- ncol(g)
  } else {
    stop("either data or graph must be provided")
  }
  if (!is.null(g)) {
    graph_edge <- c(rbind(col(g)[g == 1], row(g)[g == 1]))
    if (length(graph_edge)) graph_temp <- graph(graph_edge)
    else graph_temp <- graph.empty()
    v_order <- topological.sort(graph_temp)    
    ancestrals <- list()
    for (v in v_order) {      
      parents <- which(g[v, ] == 1)    
      ancestrals[[v]] <- union(parents, unlist(ancestrals[parents]))
      ancestrals[[v]] <- sort.int(ancestrals[[v]])      
    }    
    out$ancestrals <- ancestrals
  }
  if (is.null(g) && !is.null(model)) {
    stop("a model parameterization cannot be provided without providing the causal graph")
  }  
  if (x < 1 || x > num_v || y < 0 || y > num_v) {
    stop("treatment and/or outcome are not indexed properly")
  }
  if (!is.null(latent_idx) && (min(latent_idx) < 0 || max(latent_idx) > num_v)) {
    stop("latent variables are not indexed properly")
  }
  if (any(latent_idx) == x || any(latent_idx) == y) {
    stop("treatment and/or outcome cannot be latent variables")
  }
  
  if (num_v <= num_v_max && !is.null(model)) {
    
    probs <- array(1, rep(2, num_v))    
    for (v in v_order) {
      parents <- which(g[v, ] == 1)    
      cpt <- model[[v]]
      np <- length(parents)
      if (np > 0) {
        bef <- sum(parents < v)
        perm <- c(seq(from = np, by = -1, length = bef), np + 1, rev(seq_len(np - bef)))
        cpt <- aperm(array(c(1 - cpt, cpt), rep(2, np + 1)), perm)
      } else {
        cpt <- c(1 - cpt, cpt)
      }      
      patt <- patternRepeat(cpt, c(v, parents), rep(2, num_v))
      probs <- probs * patt
    }
    
    out$probs <- probs
    
  }  
  
  varnames <- colnames(dat)
  if (length(varnames) == 0) {
    varnames <- list()
    for (i in seq_len(num_v)) {
      varnames[[i]] <- paste("X", i, sep = "")
    }
  }
  out$varnames <- varnames
  
  class(out) <- "cfx"
  return(out)
}

#' @title Prints a CausalFX Problem Instance
#' 
#' @description 
#' Prints some of the information regarding a \code{\link{cfx}} object.
#'
#' @param x a \code{cfx} object.
#' @param ... other parameters, ignored.
#'
#' @details
#' The information that is displayed includes the identifiers of the treatment and outcome variables,
#' the names of all variables and, if a theoretical causal graph is part of the object specification,
#' its causal structured described in terms of the parent ids for each variable.
#' 
#' @export

print.cfx <- function(x, ...) {
  if (class(x) != "cfx") {
    stop("a CausalFX object is necessary")
  } 
  
  cat("CausalFX (cfx) object\n\n")
  cat(sprintf("Treatment: %s, Outcome: %s\n", x$varnames[[x$X_idx]], x$varnames[[x$Y_idx]]))
  if (is.null(x$data)) {
    cat("No dataset provided\n")
    nvar <- ncol(x$graph)
  } else {
    cat(sprintf("Dataset size: %d samples.\n", nrow(x$data)))
    nvar <- ncol(x$data)
  }
  cat("\n")
  
  is_latent <- rep(0, nvar)
  is_latent[x$latent_idx] <- 1
  
  if (is.null(x$graph)) {
    cat("Variables:\n\n")
    for (i in seq_len(nvar)) {
      cat(sprintf("    (%d) %s", i, x$varnames[[i]]))
      if (is_latent[i]) cat(" (LATENT)")
      cat("\n")
    }        
  } else {
    cat("Theoretical graph provided.\n")
    cat("Variables [with parents]:\n\n")
    for (i in seq_len(nvar)) {
      cat(sprintf("    (%d) %s", i, x$varnames[[i]]))
      if (is_latent[i]) cat(" (LATENT)")
      if (sum(x$graph[i,]) > 0) {
        cat(", parents: [")
        cat(which(x$graph[i,] == 1))
        cat("]\n")
      } else {
        cat(", parents: [NONE]\n")
      }
    }    
  }
  
}

