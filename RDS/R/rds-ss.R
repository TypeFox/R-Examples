

#' Gile's SS Estimates
#' 
#' This function computes the sequential sampling (SS) estimates for a
#' categorical variable or numeric variable.
#' 
#' 
#' @param rds.data An \code{rds.data.frame} that indicates recruitment patterns
#' by a pair of attributes named ``id'' and ``recruiter.id''.
#' @param outcome.variable A string giving the name of the variable in the
#' \code{rds.data} that contains a categorical or numeric variable to be
#' analyzed.
#' @param N An estimate of the number of members of the population being
#' sampled. If \code{NULL} it is read as the \code{population.size.mid} attribute of
#' the \code{rds.data} frame. If that is missing it defaults to 1000.
#' @param subset An optional criterion to subset \code{rds.data} by. It is a
#' character string giving an R expression which, when evaluated, subset the
#' data. In plain English, it can be something like \code{"seed > 0"} to
#' exclude seeds. It can be the name of a logical vector of the same length of
#' the outcome variable where TRUE means include it in the analysis. If
#' \code{NULL} then no subsetting is done.
#' @param number.ss.samples.per.iteration The number of samples to take in
#' estimating the inclusion probabilites in each iteration of the sequential
#' sampling algorithm. If \code{NULL} it is read as the
#' eponymous attribute of \code{rds.data}. If that
#' is missing it defaults to 5000.
#' @param number.ss.iterations The number of iterations of the sequential
#' sampling algorithm. If that is missing it defaults to 5.
#' @param control A list of control parameters for algorithm
#' tuning. Constructed using\cr
#' \code{\link{control.rds.estimates}}.
#' @param hajek logical; Use the standard Hajek-type estimator of Gile (2011)
#' or the standard Hortitz-Thompson. The default is TRUE.
#' @param empir.lik If true, and outcome.variable is numeric, standard errors
#' based on empirical likelihood will be given.
#' @return 
#' If \code{outcome.variable} is numeric then the Gile SS estimate of the mean is returned, otherwise a vector of proportion estimates is returned.
#' If the \code{empir.lik} is true, an object of class \code{rds.interval.estimate} is returned. This is a list with components
#' \itemize{ \item\code{estimate}: The numerical point estimate of proportion
#' of the \code{trait.variable}.  \item\code{interval}: A matrix with six
#' columns and one row per category of \code{trait.variable}: \itemize{
#' \item\code{point estimate}: The HT estimate of the population mean.
#' \item\code{95\% Lower Bound}: Lower 95\% confidence bound.  \item\code{95\%
#' Upper Bound}: Upper 95\% confidence bound.  \item\code{Design Effect}: The
#' design effect of the RDS.  \item\code{s.e.}: Standard error.  \item\code{n}:
#' Count of the number of sample values with that value of the trait.  } }
#' 
#' Otherwise, an object of class \code{rds.SS.estimate} is returned.
#' 
#' @author Krista J. Gile with help from Mark S. Handcock
#' @seealso \code{\link{RDS.I.estimates}}, \code{\link{RDS.II.estimates}}
#' @references
#' 
#' Gile, Krista J. 2011 \emph{Improved Inference for Respondent-Driven Sampling
#' Data with Application to HIV Prevalence Estimation}, \emph{Journal of the
#' American Statistical Association}, 106, 135-146.
#' 
#' Gile, Krista J., Handcock, Mark S., 2010 \emph{Respondent-driven Sampling:
#' An Assessment of Current Methodology}, \emph{Sociological Methodology}, 40,
#' 285-327.
#' 
#' Gile, Krista J., Handcock, Mark S., 2011 \emph{Network Model-Assisted
#' Inference from Respondent-Driven Sampling Data}, \emph{ArXiv Preprint}.
#' 
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
#' data(fauxmadrona)
#' RDS.SS.estimates(rds.data=fauxmadrona,outcome.variable="disease",N=1000)
#' 
#' @export RDS.SS.estimates
RDS.SS.estimates <-
		function(rds.data,outcome.variable,N=NULL,subset=NULL,
				number.ss.samples.per.iteration=500,number.ss.iterations=5,
				control=control.rds.estimates(),
				hajek=TRUE,empir.lik=TRUE){
	se <- substitute(subset)
	subset <- eval(se,rds.data,parent.frame())
	if(length(outcome.variable) == 1){
		result <- RDS.estimates.local(
				rds.data=rds.data,
				outcome.variable=outcome.variable,
				N=N,
				subset=subset,
				number.ss.samples.per.iteration=number.ss.samples.per.iteration,
				number.ss.iterations=number.ss.iterations,
				control=control,
				hajek=hajek,
				empir.lik=empir.lik,
				weight.type="Gile's SS")
	}
	else {
		result <- lapply(outcome.variable,function(g){
					RDS.estimates.local(
							rds.data=rds.data,
							outcome.variable=g,
							N=N,
							subset=subset,
							number.ss.samples.per.iteration=number.ss.samples.per.iteration,
							number.ss.iterations=number.ss.iterations,
							control=control,
							hajek=hajek,
							empir.lik=empir.lik,
							weight.type="Gile's SS")
					})
		names(result) <- outcome.variable
	}
	return(result)
}


