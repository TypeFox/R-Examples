#' Report quantiles and probabilities from the fitted probability distributions
#' 
#' Having fitted appropriate distributions to one or more expert's judgements
#' individually using the \code{\link{fitdist}} command, use this command to
#' get quantiles and probabilities from the fitted distributions
#' 
#' 
#' @param fit The output of a \code{fitdist} command.
#' @param quantiles A vector of desired quantiles for feedback. If this
#' argument is left out, the default is to use the same quantiles that were
#' elicited from the experts.
#' @param values A vector of desired probabilities; desired values of a for
#' reporting back fitted values of P(X<a). If this argument is left out, the
#' default is to use the same values provided by the experts.
#' @param dist If \code{fit} contains judgements from multiple experts,
#' \code{dist} is distribution to be used for calculating probabilities and
#' quantiles. Options are \code{"normal"}, \code{"t"}, \code{"gamma"},
#' \code{"lognormal"}, \code{"logt"}, \code{"beta"}, or \code{"best"}. The
#' default option, \code{"best"}, uses the best fitting distribution for each
#' expert.
#' @param ex If \code{fit} contains judgements from multiple experts,
#' specifying a value for \code{ex} will select a single expert for feedback.
#' Note that for a single expert, feedback is given for all suitable types of
#' distribution, but for multiple experts, feedback is given for one type of
#' distribution only.
#' @param sf The number of significant figures to be displayed in the output.
#' @return %% ~Describe the value returned %% If it is a LIST, use
#' 
#' \item{fitted.quantiles }{Fitted quantiles for each expert}
#' \item{fitted.probs }{Fitted probabilities for each expert}
#' \item{distributions}{The distribution used to calculate fitted
#' probabilities/quantiles for each expert, if feedback is given for multiple
#' experts.}
#' @author Jeremy Oakley <j.oakley@@sheffield.ac.uk>
#' @examples
#' \dontrun{
#' # Two experts
#' # Expert 1 states P(X<30)=0.25, P(X<40)=0.5, P(X<50)=0.75
#' # Expert 2 states P(X<20)=0.25, P(X<25)=0.5, P(X<35)=0.75
#' # Both experts state 0<X<100.
#' 
#' v <- matrix(c(30, 40, 50, 20, 25, 35), 3, 2)
#' p <- c(0.25, 0.5, 0.75)
#' myfit <- fitdist(vals = v, probs = p, lower = 0, upper = 100)
#' 
#' feedback(myfit)
#' 
#' # Feedback P(X<60) and the tertiles
#' feedback(myfit, values=60, quantiles=c(0.33,0.66))
#' 
#' # Compare fitted tertiles for different distributions, expert 2 only
#' feedback(myfit, quantiles=c(0.33,0.66), ex=2)
#' }
#' 
#' @export
feedback <- function(fit, quantiles =  NA, values = NA, dist= "best", ex = NA, sf = 3){
  if(nrow(fit$vals)>1 & is.na(ex)==T){
    return(feedbackgroup(fit, quantiles, values, dist, sf))
  }
  
  if(nrow(fit$vals)>1 & is.na(ex)==F){
    return(feedbacksingle(fit, quantiles, values, sf, ex))
  }
  
  if(nrow(fit$vals)==1){
    return(feedbacksingle(fit, quantiles, values, sf))
  }
}
