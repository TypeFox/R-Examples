#' Global test of clustering using difference in K functions
#' 
#' \code{kdplus.test} performs a global test of clustering for comparing cases and controls using the method of Diggle and Chetwynd (1991).  It relies on the difference in estimated K functions.
#' 
#' @param x A \code{kdenv} object from the \code{kdest} function.  
#' 
#' @return A list providing the observed test statistic (\code{kdplus}) and the estimate p-value {\code{pvalue}}.
#' @author Joshua French
#' @import spatstat
#' @importFrom stats sd
#' @seealso \code{\link{kdest}}
#' @export
#' @references Waller, L.A. and Gotway, C.A. (2005).  Applied Spatial Statistics for Public Health Data.  Hoboken, NJ: Wiley.  Diggle, Peter J., and Amanda G. Chetwynd. "Second-order analysis of spatial clustering for inhomogeneous populations." Biometrics (1991): 1155-1163.
#' @examples 
#' data(grave)
#' kdsim = kdest(grave, nsim = 9)
#' kdplus.test(kdsim)

kdplus.test = function(x)
{
  if(max(class(x) == "kdenv") < 1) stop("x must be an object from the kdenv function.")
  simfuns <- as.data.frame(attr(x[[1]], "simfuns"))
  simfuns[,1] <- x[[1]]$obs # replace r with obs kd
  sdkdhat = apply(simfuns, 1, stats::sd) # estimated variance of kdest simulations
  # turn into matrix
  sdmat = matrix(sdkdhat, nrow = nrow(simfuns), ncol = ncol(simfuns))
  # estimate KD + for simulated data
  kdplussim = colSums(simfuns/sdmat, na.rm = TRUE)
  # determine proportion of simulated KD+ and observed KD+
  # greater than KD+
  p = mean(kdplussim >= kdplussim[1])
  print(paste("The p-value for the global test is", round(p, 3)))
  return(invisible(list(kdplus = kdplussim[1], pvalue = p)))
}