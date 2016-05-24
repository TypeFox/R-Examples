#' Propose a new real hyperparameter value.
#' 
#' This function proposes a new real values hyperparameter for the information
#' sharing prior.
#' 
#' 
#' @param orig_beta Current value of the hyperparameter.
#' @param proposal_range Range for the new value.
#' @param limit Hard limit on the range.
#' @return Returns a new uniformly random value within \code{proposal_range} of
#' \code{orig_beta} and limited by \code{limit}.
#' @author Frank Dondelinger
#' @seealso \code{\link{ProposeDiscrete}}
#' @examples
#' 
#' # Previous parameter value
#' param = runif(1, 0, 1)
#' 
#' # Propose new value within range [0, 1], with proposal width 0.1
#' new.param = proposeContinuous(param, 0.1, 1)
#' 
#' @export proposeContinuous
proposeContinuous <-
function(orig_beta, proposal_range, limit=30) {

  a = orig_beta - proposal_range;
  b = orig_beta + proposal_range;

  new_beta = a + (b-a)*runif(1, 0, 1);
 
  if (new_beta > limit)
    new_beta = (2*limit) - new_beta;
 
  return(abs(new_beta))

}

