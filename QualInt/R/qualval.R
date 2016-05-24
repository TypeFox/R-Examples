#' Test for qualitative interactions from estimation
#' 
#' Test for qualitative interactions between treatment effects and patient subgroups
#' from the estimated treatment effect and its standard error directly.
#' Output all the results related with qualitative interaction tests as 
#' a "qualint" object, just like the "qualint" function. 
#' Two common tests for qualitative interactions are included: IBGA and LRT, among
#' which IBGA is the default. 
#' Users need to input the estiamted treatment effect and its standard erro themselves,
#' therefore could accommodate any types of responses.
#' 
#' @param effect treatment effects. A numeric vector as the estimated (for pvalue) or 
#' designed (for power) treatment effects for different patient subgroups.
#' @param se standard error of estimated treatmen effects. A numeric vector as the
#' standard error for the estimated treatment effects.
#' @param test testing method. Choose either \code{"IBGA"} (interval based graphical
#' approach) or \code{"LRT"} (Gail Simon likelihood ratio test).
#' @param alpha significance level. The type I error for qualitative 
#' interaction tesing. The default is 0.05.
#' @param plotout whether output the plot or not for \code{test = "IBGA"}.
#' There is no plot output for \code{test = "LRT"}.
#' 
#' @details
#' This function is a more generalized version of \code{qualint} in the sense that
#' it could be used for any types of responses. However, comepared to \code{qualint},
#' the user needs to input the estimated (for pvalue) or designed (for power ) 
#' treatment effects and its standard error by themselves to use this function. It
#' gives more freedom and allows users to choose the method they prefer before
#' testing for qualitative interactions.
#' 
#' @return An object with S3 class "qualint".
#' \item{call}{the call that produces this object.}
#' \item{n}{the sample size for each treatment in each subgroup.}
#' \item{type}{response type.}
#' \item{alpha}{significance level for the test.}
#' \item{treatment}{treatment factors.}
#' \item{reference}{reference treatment used for the comparison.}
#' \item{nsbp}{the number of patient subgroups.}
#' \item{subgroup}{subgroup factors.}
#' \item{scale}{the scale type for treatment effects (see above). }
#' \item{effect}{estimated treatment effects.}
#' \item{se}{standard error of treatment effects estimators.}
#' \item{LowerCI}{the lower limit of the confidence interval.}
#' \item{UpperCI}{the upper limit of the confidence interval.}
#' \item{test}{testing method used here, either "IBGA" or "LRT".}
#' \item{index}{the testing index used only for \code{test = "IBGA"}.}
#' \item{cvalue}{the critical value used only for \code{test = "LRT"}.}
#' \item{LowerTI}{the lower limit of the testing interval used when \code{test = "IBGA"}.}
#' \item{UpperTI}{the upper limit of the testing interval used when \code{test = "IBGA"}.}
#' \item{pvalue}{the pvalue for qualitative interactions.}
#' \item{power}{the power based on the observed data.}
#' \item{nobs}{the number of subjects.}
#' \item{missing}{the indexes of subjects with missing values.}
#' 
#' @author Lixi Yu, Eun-Young Suh, Guohua (James) Pan \cr
#' Maintainer: Lixi Yu \email{lixi-yu@@uiowa.edu}
#' 
#' @references Gail and Simon (1985), Testing for qualitative interactions between
#' treatment effects and patient subsets, Biometrics, 41, 361-372.
#' @references Pan and Wolfe (1993), Tests for generalized problems of detecting
#' qualitative interaction, Technical Report No. 526, Department of Statistics,
#' The Ohio State University.
#' @references Pan and Wolfe (1997), Test for qualitative interaction of clinical
#' significance, Statistics in Medicine, 16, 1645-1652. 
#' 
#' @seealso \code{\link{print.qualint}}, \code{\link{coef.qualint}}, 
#' \code{\link{plot.qualint}}, \code{\link{qualint}}
#' 
#' @examples
#' 
#' test9 <- qualval(effect = c(1.0, 0.5, -2.0), se = c(0.86, 0.64, 0.32))
#' print(test9)
#' plot(test9)
#' 
#' @export

qualval <- function(effect, se, test = c("IBGA", "LRT"), 
                    alpha = 0.05, plotout = FALSE) {
  
  test <- match.arg(test)
  
  alpha <- as.vector(alpha) 
  if(length(alpha) != 1) {
    warning("alpha has length > 1 and only the first element will be used")
    alpha <- alpha[1]
  }   
  if(is.null(alpha)) stop("argument alpha can not be NULL")
  alpha <- as.numeric(alpha)
  if(is.na(alpha) | (alpha >= 1) | (alpha <= 0))
    stop("alpha should be a numeric value between 0 and 1")
  
  effect <- as.vector(effect)
  se <- as.vector(se)
  if(!(is.numeric(effect) & is.numeric(se)))
    stop("arguments effect and se should both be numeric")
  if(length(effect) != length(se))
    stop("the length of argument effect and se should be the same")
  
  nsbp <- length(effect)
  betaout <- matrix(c(effect, se), nsbp, 2)
  testout <- switch(test, 
                    IBGA = ibgaout(betaout, nsbp, alpha),
                    LRT = lrtout(betaout, nsbp, alpha))
  
  this_call = match.call()
  CI_index <- qnorm(1 - alpha / 2)
  obj <- switch(test,
                IBGA = list(call = this_call, n = NULL, type = "unknown", 
                            alpha = alpha, treatment = "1", reference = "0", 
                            nsbp = nsbp, subgroup = as.character(1 : nsbp), 
                            scale = "unknown", effect = betaout[, 1], se = betaout[, 2],
                            LowerCI = betaout[, 1] - CI_index * betaout[, 2],
                            UpperCI = betaout[, 1] + CI_index * betaout[, 2],
                            test = test, index = testout$index, 
                            LowerTI = betaout[, 1] - testout$index * betaout[, 2],
                            UpperTI = betaout[, 1] + testout$index * betaout[, 2],
                            pvalue = testout$pvalue, power = testout$power, nobs = NULL,
                            missing = NULL),
                LRT = list(call = this_call, n = NULL, type = "unknown", 
                           alpha = alpha, treatment = "1", reference = "0", 
                           nsbp = nsbp, subgroup = as.character(1 : nsbp), 
                           scale = "unknown", effect = betaout[, 1], se = betaout[, 2],
                           LowerCI = betaout[, 1] - CI_index * betaout[, 2],
                           UpperCI = betaout[, 1] + CI_index * betaout[, 2],
                           test = test, cvalue = testout$index, 
                           pvalue = testout$pvalue, power = testout$power, nobs = NULL, 
                           missing = NULL))
  
  class(obj) <- "qualint"
  
  plotout <- as.vector(plotout)
  if(length(plotout) != 1) {
    warning("plot has length > 1 and only the first element will be used")
    plotout <- plotout[1]
  } 
  if(!is.logical(plotout)){
    warning("plot should be a logical value and set to FALSE")
    plotout <- FALSE
  }
  if(plotout & test == "LRT")
    warning("there is no plot output for LRT")
  
  if(plotout & test == "IBGA") plot(obj)
  obj
  
}