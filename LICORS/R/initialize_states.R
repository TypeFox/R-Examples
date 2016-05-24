#' @title State initialization for iterative algorithms (randomly or variants of kmeans)
#'
#' @description 
#' Initializes the state/cluster assignment either uniformly at random from
#' \eqn{K} classes, or using initial \emph{kmeans++} (\code{\link{kmeanspp}}) 
#' clustering (in several variations on PLCs and/or FLCs).
#' 
#' @param num.states number of states 
#' @param num.samples number of samples.
#' @param method how to choose the labels: either uniformly at random from 
#' \eqn{\lbrace 1, \ldots, K \rbrace} or using K-means on PLCs and FLCs or 
#' a combination.  Default: \code{method = "random"}.  Other options are
#' \code{c("KmeansPLC","KmeansFLC","KmeansPLCFLC","KmeansFLCPLC")}
#' @param LCs (optional) a list of \code{PLC} (\eqn{N \times n_p} array) and 
#' \code{FLC} (\eqn{N \times n_f} array)
#' @keywords datagen distribution multivariate
#' @export
#' @examples
#' x1 = rnorm(1000)
#' x2 = rnorm(200, mean = 2)
#' yy = c(x1, x2)
#' ss = initialize_states(num.states = 2, num.samples = length(yy), method = "KmeansFLC", 
#'                        LCs = list(FLCs = yy))
#' plot(yy, col = ss, pch = 19)
#' points(x1, col = "blue")
#' 

initialize_states <- function(num.states = NULL, num.samples = NULL, 
                              method = c("random", "KmeansPLC", "KmeansFLC", 
                                         "KmeansPLCFLC", "KmeansFLCPLC"),
                              LCs = list(PLC = NULL, FLC = NULL)) {
  if (is.null(num.states)) {
    stop("You must provide the number of clusters.")
  }
  
  if (is.null(num.samples)) {
    if (is.null(unlist(LCs))) {
      stop("You must either provide the total number of samples or data 'LCs'.")
    } else {
      num.samples <- nrow(LCs$PLC)
    }
  }  
  method <- match.arg(method)

  if (method == "random") {
    states <- sample.int(n = num.states, size = num.samples, replace = TRUE)
    # make sure that every state at least appears once
    states[sample.int(num.samples, num.states, replace = FALSE)] <- 
      seq_len(num.states)
  } else if (method == "KmeansPLC") {
    states <- kmeanspp(LCs$PLC, num.states, iter.max = 100, nstart = 10)$cluster
  } else if (method == "KmeansFLC") {
    states <- kmeanspp(LCs$FLC, num.states, iter.max = 100, nstart = 10)$cluster
  } else if (any(method == c("KmeansPLCFLC", "KmeansFLCPLC"))) {
    # do two stage clustering 1) on the PLC (FLC) space, and then conditional
    # on the first stage cluster 2) do a clustering in each cluster but using
    # the FLCs (PLCs)
    first.stage.num.states <- ceiling(num.states^(1/3))
    second.stage.num.states <- floor(num.states/first.stage.num.states)
    second.stage.last.run.num.states <- 
      num.states - first.stage.num.states * second.stage.num.states + second.stage.num.states
    
    if (method == "KmeansFLCPLC") {
      first.stage.data <- rbind(LCs$FLC)
      second.stage.data <- rbind(LCs$PLC)
    } else {
      first.stage.data <- rbind(LCs$PLC)
      second.stage.data <- rbind(LCs$FLC)
    }
        
    first.stage.labels <- kmeanspp(first.stage.data, first.stage.num.states, 
                                   iter.max = 100, nstart = 10)$cluster

    second.stage.labels <- rep(NA, num.samples)
    for (ll in seq_len(first.stage.num.states - 1)) {
      second.stage.labels[first.stage.labels == ll] <- (ll - 1) * first.stage.num.states + 
        kmeanspp(second.stage.data[first.stage.labels == ll,], 
                 second.stage.num.states, 
                 iter.max = 100, nstart = 10)$cluster
    }
    second.stage.labels[first.stage.labels == first.stage.num.states] = (first.stage.num.states - 1) * first.stage.num.states + 
      kmeanspp(second.stage.data[first.stage.labels == first.stage.num.states,], 
               second.stage.last.run.num.states, 
               iter.max = 100, nstart = 10)$cluster
    
    states <- second.stage.labels
  }
  invisible(states)
}

