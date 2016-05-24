#' Propose a new discrete value.
#' 
#' This function proposes a new discrete parameter, based on the previous
#' value, within the given proposal range, making sure that the maximum range
#' is not exceeded.
#' 
#' 
#' @param params.old Old parameter value (an integer).
#' @param proposal.range Range for new proposal (an integer).
#' @param max.range Maximum value for new proposal (an integer).
#' @return Returns the new proposed parameter, which will be an integer in the
#' range [0, \code{max.range}], and within at most \code{proposal.range} of
#' \code{params.old}.
#' @author Frank Dondelinger
#' @seealso \code{\link{proposeContinuous}}
#' @examples
#' 
#' # Previous parameter value
#' param = rpois(1, 5)
#' 
#' # Propose new value within range [0, 10], with proposal width 2
#' new.param = ProposeDiscrete(param, 2, 10)
#' 
#' @export ProposeDiscrete
ProposeDiscrete <-
function(params.old, proposal.range, max.range) {
  # Propose a new discrete parameter, based on the previous value, within the
  # given proposal range, making sure that the maximum range is not exceeded. 
  #
  # Args:
  #  params.old:     Old value for the parameter
  #  proposal.range: Range for the proposal
  #  max.range:      Total range for the parameter
  #
  # Returns:
  #  New value for the parameter
  
  extent = (params.old - proposal.range):(params.old + proposal.range)
  
  pick = sample(1:length(extent), 1)
  
  params.new = extent[pick]
  
  if(params.new < 1) { 
    params.new = max.range + params.new
  } else if(params.new > max.range) {
    params.new = params.new - max.range
  }
  
  return(params.new)
}

