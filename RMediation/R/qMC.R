#' Quantile for the Monte Carlo Sampling Distribution of a nonlinear function of coefficients estimates
#'
#' This function returns a quantile corresponding to the probability \code{p}.
#'
#'@param p probability.
#' @param mu a \link{vector} of means (e.g., coefficient estimates) for the normal random variables. A user can assign a name to each mean value, e.g., \code{mu=c(b1=.1,b2=3)}; otherwise, the coefficient names are assigned automatically as follows: \code{b1,b2,...}.
#' @param Sigma either a covariance matrix or a \link{vector} that stacks all the columns of the lower triangle variance--covariance matrix one underneath the other.
#' @param quant quantity of interest, which is a nonlinear/linear function of the model parameters. Argument \code{quant} is a \link{formula} that \strong{must} start with the symbol "tilde" (\code{~}): e.g., \code{~b1*b2*b3*b4}. The names of coefficients must conform to the names provided in the argument \code{mu} or to the default names, i.e., \code{b1,b2,...}.
#' @param n.mc Monte Carlo sample size. The default sample size is 1e+6.
#' @param ... additional arguments.
#' @return scalar quantile value.
#' @keywords regression distribution
#' @export qMC
#' @examples
#' qMC(.05,mu=c(b1=1,b2=.7,b3=.6, b4= .45), Sigma=c(.05,0,0,0,.05,0,0,.03,0,.03), 
#' quant=~b1*b2*b3*b4)
#' @author Davood Tofighi \email{dtofighi@@psych.gatech.edu} and David P. MacKinnon \email{davidpm@@asu.edu}
#' @references  Tofighi, D. and MacKinnon, D. P. (2011). RMediation: An R package for mediation analysis confidence intervals. \emph{Behavior Research Methods}, \bold{43}, 692--700. doi:10.3758/s13428-011-0076-x
#' @seealso \code{\link{medci}} \code{\link{RMediation-package}}


qMC <- function(p, mu, Sigma, quant, n.mc = 1e+06,...){
  if(missing(p)) stop(paste("argument",sQuote("p"), "must be specified"))
  if(missing(mu)) stop(paste("argument",sQuote("p"), "must be specified"))
  if(missing(Sigma)) stop(paste("argument",sQuote("Sigma"), "must be specified"))
  if(missing(quant)) stop(paste("argument",sQuote("quant"), "must be specified"))
  if(is.null(p)) stop(paste("argument",sQuote("mu"), "cannot be a NULL value"))
  
  if(is.null(mu)) stop(paste("argument",sQuote("mu"), "cannot be a NULL value"))
  if(is.null(Sigma)) stop(paste("argument",sQuote("Sigma"), "cannot be a NULL value"))
  if(is.null(quant)) stop(paste("argument",sQuote("NULL"), "cannot be a NULL value"))
   if(!is.matrix(Sigma)){
    if(length(mu)!= (sqrt(1 + 8 * length(Sigma)) - 1)/2) stop(paste("Please check the length of", sQuote("Sigma"),"and",sQuote("mu"),". If the length(dimension) of the", sQuote("mu"),"vector (",length(mu),") is correct, the stacked lower triangle matrix", sQuote("Sigma"), "must have ",((2*length(mu)+1)^2-1)/8, "elements, instead of", length(Sigma)) )
    Sigma <- lav_matrix_vech_reverse(Sigma) #converts to a symmetric matrix
  }
  if(is.null(names(mu)) ) names(mu) <- paste("b",1:length(mu), sep="") # if mu names is NULL
  
  if(!all(all.vars(quant) %in% names(mu))) stop(paste("The parameters names in formula", sQuote("quant"), "must match the parameters names provided in", sQuote("mu"),"."))
  
  if(length(mu)*n.mc > .Machine$integer.max){
    n.mc <- 1e6
    warning(paste("n.mc is too large. It is reset to", n.mc))
  }
  if(p<=0 |p>=1) stop('Enter a correct p value, 0<p<1') 
  
  quant <- parse(text=sub("~","",quant))
  df <- data.frame(mvrnorm(n.mc,mu,Sigma))
  colnames(df) <-names(mu)
  quant.vec <- eval(quant,df)
  return(quantile(quant.vec,p) )
 }

