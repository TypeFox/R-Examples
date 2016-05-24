#' @title Estimate conditional/marginal state probabilities
#'
#' @description 
#' Estimates \eqn{P(S = s_k; \mathbf{W})}, \eqn{k = 1, \ldots, K}, the probability of being 
#' in state \eqn{s_k} using the weight matrix \eqn{\mathbf{W}}.
#' 
#' These probabilites can be marginal (\eqn{P(S = s_k; \mathbf{W})}) or 
#' conditional (\eqn{P(S = s_k \mid \ell^{-}, \ell^{+}; \mathbf{W})}), depending on the 
#' provided information (\code{pdfs$PLC} and \code{pdfs$FLC}).  
#' \itemize{
#'   \item If both are \code{NULL} then \code{estimate_state_probs} returns a 
#'         vector of length \eqn{K} with marginal probabilities.
#'   \item If either of them is not \code{NULL} then it returns an 
#'         \eqn{N \times K} matrix, where row \eqn{i} is the 
#'         probability mass function of PLC \eqn{i} being in state 
#'         \eqn{s_k}, \eqn{k = 1, \ldots, K}.
#' }
#'  
#' @param weight.matrix \eqn{N \times K} weight matrix
#' @param states vector of length \eqn{N} with entry \eqn{i} being the label 
#' \eqn{k = 1, \ldots, K} of PLC \eqn{i}
#' @param pdfs a list with estimated pdfs for PLC and/or FLC evaluated at each PLC, \eqn{i=1, \ldots, N}
#' and/or FLC, \eqn{i=1, \ldots, N}
#' @param num.states number of states in total. If \code{NULL} (default) then
#' it sets it to \code{max(states)} or \code{ncol(weight.matrix)} - 
#' depending on which one is provided.
#' @return
#' A vector of length \eqn{K} or a \eqn{N \times K} matrix.
#' @keywords nonparametric multivariate distribution
#' @export
#' @examples
#' WW = matrix(runif(10000), ncol = 10)
#' WW = normalize(WW)
#' estimate_state_probs(WW)

estimate_state_probs <- function(weight.matrix = NULL, 
                                 states = NULL,
                                 pdfs = list(FLC = NULL, PLC = NULL),
                                 num.states = NULL) {
  
  if (is.null(states) & is.null(weight.matrix)) {
    stop("You must provide either a vector of state 
          labels or the weight matrix.")
  }
  if (is.null(num.states)) {
    if (is.null(states)) {
      num.states <- ncol(weight.matrix)
    } else {
      num.states <- max(states)
    }
  }
  
  if (is.null(pdfs$PLC) && is.null(pdfs$FLC)) {
    ##############################
    ### marginal probabilities ###
    ##############################    
    if (is.null(weight.matrix)) {
      state.probs <- 
        c((table(c(states, seq_len(num.states))) - 1)) / length(states)
    } else {
      state.probs <- colMeans(weight.matrix)
    } 
  } else { 
    #################################
    ### conditional probabilities ###
    #################################
    num.LCs <- ifelse(is.null(nrow(pdfs$FLC)), nrow(pdfs$PLC), nrow(pdfs$FLC))
    # initialize pi_theta
    state.probs.given.LCs <-
      Matrix(0, ncol = num.states, nrow = num.LCs, sparse = TRUE)

    # if PLCs are not provided, initialize weights with 0/1 assignment
    # based on state
    if (is.null(pdfs$PLC)) {
      state.probs.given.LCs[cbind(seq_len(nrow(state.probs.given.LCs)),
                                  states)] <- 1
    } else {
      # estimate marginal probabilities ('prior' probabilities for each state)
      marginal.state.probs <- estimate_state_probs(weight.matrix = weight.matrix, 
                                                   states = states,
                                                   num.states = num.states)
      state.probs.given.LCs <- pdfs$PLC %*% Diagonal(x = marginal.state.probs)
        #sweep(pdfs$PLC, 2, marginal.state.probs, "*")  makes sparse matrix dense
    } 
    
    # condition on FLCs if provided
    if (!is.null(pdfs$FLC)) {
      state.probs.given.LCs <- state.probs.given.LCs * pdfs$FLC
    }
  }

  if (is.null(pdfs$PLC) && is.null(pdfs$FLC)) {
    return(state.probs)
  } else {
    invisible(normalize(state.probs.given.LCs))
  }
} 
