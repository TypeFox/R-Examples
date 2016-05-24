#' Johnson Transformation for Normality
#' 
#' \code{jtrans} transforms a continuous univariate vector to a random vector 
#' from standard normal distribution.
#' 
#' \code{jtrans} fits data to a set of distributions from Johnson family. A 
#' normality test is used to find the best fit by choosing the fit with maximum 
#' p.value under that given test. It returns the transformed data, the 
#' corresponding type of Johnson curve and parameter estimations.
#' 
#' Since the default Shapiro-Wilk test can only accept sample size between 3 and
#' 5000, one should specify another normality test in the test parameter, 
#' generally the \code{\link[nortest]{ad.test}} in the \pkg{nortest} package is 
#' recommended.
#' 
#' Sometimes, this algorithm may return poor fits. The most extreme case is that
#' all the transformed data have smaller p.values than the original data's. In 
#' such cases, the \code{exclude_original} flag should be set to FALSE, so 
#' \code{jtrans} will return the original data as the transformed data.
#' 
#' @param x the non-normal numerical data.
#' @param test the normality test used to select fits, defaults to 
#'   \code{\link{shapiro.test}}
#' @param exclude_original whether the original data should be excluded when 
#'   comparing fits.
#' @param z_lim two values vector defining the range of the z values, defaults 
#'   to 0.25 to 1.25, which is recommended by Mandraccia, Halverson and Chou 
#'   (1996).
#' @param z_length the length of the z vector, default to 101. The number of 
#'   different fits estiamted in the algorithm. Set larger z.length value if you
#'   want extra precision.
#'   
#' @export
#' 
#' @return A list with two classes: the first one is the type of transformation
#'   used, the same as the \code{type} component, could be "sb", "su" or "sl";
#'   The second one is "jtrans". The list containsthe following components: 
#'   \item{original}{original data.} \item{transformed}{transformed data.} 
#'   \item{type}{type of transformation selected.} \item{test}{normality test 
#'   used to select transformations.} \item{z}{selected z value among 101 values
#'   from 0.25 to 1.25.} \item{eta, gamma, lambda, epsilon}{transformation 
#'   parameters.} \item{p.value}{the maximum p.value returned by test}
#'   
#' @references Chou, Y. M., Polansky, A. M., & Mason, R. L. (1998). Transforming
#'   non-normal data to normality in statistical process control. Journal of 
#'   Quality Technology, 30(2), 133-141.
#'   
#' @examples
#' # generate 100 non-normal data and transform it.
#' x <- rexp(50, .2)
#' jt <- jtrans(x)
#' jt
#' 
#' plot(density(x))
#' plot(density(jt$transformed))
#' qqnorm(jt$transformed)
#' qqline(jt$transformed)
#' 
#' \dontrun{
#' # Using another normality test
#' require(nortest)
#' jtrans(x, test = "ad.test")
#' }

jtrans <- function(x, test="shapiro.test", exclude_original=TRUE, 
                   z_lim = c(.25, 1.25),
                   z_length = 101) {
  # set trans to be the original data, collect original parameters
  x <- as.numeric(x)
  trans <- x
  type <- "original"
  params <- data.frame(eta     = NA, 
                       gamma   = NA, 
                       lambda  = NA, 
                       epsilon = NA,
                       z       = NA,
                       p.value = do.call(test, list(trans))$p.value)
  
  
  # test for transformed dataset corresponding to 101 z value
  for (z in seq(from=z_lim[1], to=z_lim[2], length.out = z_length)) {
    # calculate quantiles
    q <- qtls(x, z)
    
    # fit sl distribution for every z
    res <- fit_sl(x, q)
    if (!is.na(res[1])) {
      # fit succeeds
      type <- c(type, "sl")
      params <- rbind(params, c(res$params,
                                do.call(test, list(res$trans))$p.value))
      trans <- rbind(trans, res$trans)
    }
    
    # fit either SB or SL for every z accorrding to QR
    if (q$QR <= 1) { 
      res <- fit_sb(x, q)
      if (!is.na(res[1])) {
        type <- c(type, "sb")
        params <- rbind(params, c(res$params,
                                  do.call(test, list(res$trans))$p.value))
        trans <- rbind(trans, res$trans)
      }
    } else {
      res <- fit_su(x, q)
      if (!is.na(res[1])) {
        type <- c(type, "su")
        params <- rbind(params, c(res$params, 
                                  do.call(test, list(res$trans))$p.value))
        trans <- rbind(trans, res$trans)
      }
    }
  }
  
  if(exclude_original) {
    type <- type[-1]
    trans <- trans[-1, ]
    params <- params[-1, ]
  }
  
  max <- which.max(params$p.value)
  # collect max p.value
  RVAL <- structure(c(list(original    = x,
                           transformed = trans[max, ], 
                           type        = type[max],
                           test        = test),
                      params[max, ]), 
                    class = c(type[max], "jtrans"))
  return(RVAL)
}