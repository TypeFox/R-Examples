

#' Compute RDS-I Estimates
#' 
#' This function computes the RDS-I type estimates for a categorical variable.
#' It is also referred to as the Salganik-Heckathorn estimator.
#' 
#' 
#' @aliases RDS.I.estimates
#' @param rds.data An \code{rds.data.frame} that indicates recruitment patterns
#' by a pair of attributes named ``id'' and ``recruiter.id''.
#' @param outcome.variable A string giving the name of the variable in the
#' \code{rds.data} that contains a categorical variable to be analyzed.
#' @param subset An expression defining a subset of rds.data.
#' @param smoothed Logical, if TRUE then the ``data smoothed'' version of RDS-I is used,
#' where it is assumed that the observed Markov process is reversible.
#' @param empir.lik Should confidence intervals be estimated using 
#' empirical likelihood.
#' @param N Population size to be used to calculate the empirical likelihood interval. If NULL, this value is
#' taken to be the population.size.mid attribute of the data and if that is not set, no finite population
#' correction is used.
#' @return If the \code{empir.lik} is true, an object of class
#' \code{rds.interval.estimate} is returned. This is a list with components
#' \itemize{ \item\code{estimate}: The numerical point estimate of proportion
#' of the \code{trait.variable}.  \item\code{interval}: A matrix with six
#' columns and one row per category of \code{trait.variable}: \itemize{
#' \item\code{point estimate}: The HT estimate of the population mean.
#' \item\code{95\% Lower Bound}: Lower 95\% confidence bound.  \item\code{95\%
#' Upper Bound}: Upper 95\% confidence bound.  \item\code{Design Effect}: The
#' design effect of the RDS.  \item\code{s.e.}: Standard error.  \item\code{n}:
#' Count of the number of sample values with that value of the trait.  } }
#' 
#' Otherwise an object of class \code{rds.I.estimate} object is returned.
#' 
#' @author Mark S. Handcock and W. Whipple Neely
#' @seealso \code{\link{RDS.II.estimates}}, \code{\link{RDS.SS.estimates}}
#' @references Gile, Krista J., Handcock, Mark S., 2010,
#' \emph{Respondent-driven Sampling: An Assessment of Current Methodology}.
#' Sociological Methodology 40, 285-327.
#' 
#' Neely, W. W., 2009. \emph{Bayesian methods for data from respondent driven
#' sampling}. Dissertation in-progress, Department of Statistics, University of
#' Wisconsin, Madison.
#' 
#' Salganik, M., Heckathorn, D. D., 2004. \emph{Sampling and estimation in
#' hidden populations using respondent-driven sampling}. Sociological
#' Methodology 34, 193-239.
#' 
#' Volz, E., Heckathorn, D., 2008. \emph{Probability based estimation theory
#' for Respondent Driven Sampling}. The Journal of Official Statistics 24 (1),
#' 79-97.
#' @keywords survey manip
#' @examples
#' 
#' data(faux)
#' RDS.I.estimates(rds.data=faux,outcome.variable='X')
#' RDS.I.estimates(rds.data=faux,outcome.variable='X',smoothed=TRUE)
#' @export 
RDS.I.estimates <- function(rds.data,outcome.variable, N=NULL,subset=NULL, smoothed=FALSE,
		empir.lik=TRUE){
	se <- substitute(subset)
	subset <- eval(se,rds.data,parent.frame())
	weight.type <- if(smoothed) "RDS-I (DS)" else "RDS-I"
	if(length(outcome.variable) == 1){
		result <- RDS.estimates.local(rds.data,outcome.variable,
				subset=subset, empir.lik=empir.lik,weight.type=weight.type, N=N)
	}
	else {
		result <- lapply(outcome.variable,function(g){
					RDS.estimates.local(rds.data,g,subset=subset,
							 empir.lik=empir.lik,weight.type=weight.type, N=N)
		})
		names(result) <- outcome.variable
		
	}
	return(result)
}
