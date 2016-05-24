#' @title Tools to Support the Sheffield Elicitation Framework
#' 
#' @description Implements various methods for eliciting a probability distribution for a
#' single parameter from an expert or a group of experts. The expert provides a
#' small number of probability or quantile judgements, corresponding to points
#' on his or her cumulative distribution function. A range of parametric
#' distributions can then be fitted and displayed, with feedback provided in
#' the form of additional quantiles. A graphical interface for the roulette
#' elicitation method is also provided. For multiple experts, a weighted linear
#' pool can be calculated.
#' 
#' \tabular{ll}{ Package: \tab SHELF\cr Type: \tab Package\cr Version: \tab
#' 1.0.1.9002\cr Date: \tab 2016-01-11\cr License: \tab GPL-2 | GPL-3\cr }
#' 
#' @name SHELF-package
#' @aliases SHELF-package SHELF
#' @docType package
#' @author Jeremy Oakley <j.oakley@@sheffield.ac.uk>
#' @references \href{http://www.tonyohagan.co.uk/shelf/}{The SHELF homepage}
#' @examples
#' \dontrun{
#' ## 1) Elicit judgements from two experts individually 
#' # Expert 1 states P(X<30)=0.25, P(X<40)=0.5, P(X<50)=0.75
#' # Expert 2 states P(X<20)=0.25, P(X<25)=0.5, P(X<35)=0.75
#' # Both experts state 0<X<100.
#' 
#' ## 2) Fit distributions to each expert's judgements
#' v <- matrix(c(30, 40, 50, 20, 25, 35), 3, 2)
#' p <- c(0.25, 0.5, 0.75)
#' myfit <- fitdist(vals = v, probs = p, lower = 0, upper = 100)
#' 
#' ## 3) Plot the fitted distributions, including a linear pool
#' plotfit(myfit, lp=T)
#' 
#' ## 4) Now elicit a single 'consensus' distribution from the two experts
#' # Suppose they agree P(X<25)=0.25, P(X<30)=0.5, P(X<40)=0.75
#' v <-c(25, 30, 40)
#' p <-c(0.25, 0.5, 0.75)
#' myfit <- fitdist(vals = v, probs = p, lower = 0, upper = 100)
#' 
#' ## 5) Plot the fitted density, and report some feedback, such as the 
#' # fitted 5th and 95th percentiles
#' plotfit(myfit, ql=0.05, qu=0.95)
#' feedback(myfit, quantiles=c(0.05,0.95))
#' 
#' ## Can also use interactive plotting
#' v <- matrix(c(30, 40, 50, 20, 25, 35), 3, 2)
#' p <- c(0.25, 0.5, 0.75)
#' myfit <- fitdist(vals = v, probs = p, lower = 0, upper = 100)
#' # plot each distribution
#' plotfit(myfit, int = TRUE)
#' 
#' ## plot the distribution for one expert only
#' plotfit(myfit, int= TRUE, ex=1)
#' 
#' ## Enter judgements in interactive mode
#' elicit()
#' 
#' ## Enter judgements using the roulette method
#' roulette(Lo=0, Up=100, nbins=10, gridheight=10)
#' }
{
}

