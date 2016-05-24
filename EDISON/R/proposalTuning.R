#' Tune the proposal width for betas.
#' 
#' This function adjusts the proposal width for the beta hyperparameter(s) of
#' the exponential information sharing prior, so that the acceptance rate is
#' close to 0.25.
#' 
#' 
#' @param acceptRate Current acceptance rate.
#' @param hyper.proposals Current proposal width.
#' @return Returns the new proposal width.
#' @author Frank Dondelinger
#' @export proposalTuning
proposalTuning <-
function(acceptRate, hyper.proposals) {
  if(is.nan(acceptRate[2])) return(hyper.proposals)
  
  decrement = 1 - runif(1, 0, 0.5)
  increment = 1 + runif(1, 0, 0.5)
  
  if(acceptRate[2] < 0.2) {
    hyper.proposals = hyper.proposals*decrement
  } else if(acceptRate[2] > 0.3) {
    hyper.proposals = hyper.proposals*increment
  }

  return(hyper.proposals)
}

