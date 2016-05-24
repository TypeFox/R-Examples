#' Global test of clustering using log ratio of spatial densities
#' 
#' \code{logrr.test} performs a global test of clustering for comparing cases and controls using the log ratio of spatial densities the method of Kelsall and Diggle (1995).
#' 
#' @param x An \code{renv} object from the \code{logrr.env} function.
#' 
#' @return A list providing the observed test statistic (\code{islogrr}) and the estimated p-value {\code{pvalue}}.
#' @author Joshua French
#' @import spatstat
#' @export
#' @references Waller, L.A. and Gotway, C.A. (2005).  Applied Spatial Statistics for Public Health Data.  Hoboken, NJ: Wiley.  Kelsall, Julia E., and Peter J. Diggle. "Non-parametric estimation of spatial variation in relative risk." Statistics in Medicine 14.21-22 (1995): 2335-2342.
#' @examples 
#' data(grave)
#' renv = logrr(grave, nsim = 9)
#' logrr.test(renv)

logrr.test = function(x)
{
  if(max(class(x) == "logrrenv") < 1) stop("x must be an object from the logrr function.")
  
  win = x$window
  dim3 = dim(x$simr)[3]
  stats = numeric(dim3)
  for(i in 1:dim3)
  {
   x$v = x$simr[,,i]^2
   stats[i] = integral.im(x)
  }
  islogrr = stats[1]
  p = mean(stats >= islogrr)
  print(paste("The p-value for the global test is", p))
  return(invisible(list(islogrr = islogrr, pvalue = p)))
}