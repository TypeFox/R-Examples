#' Test for qualitative interactions from complete data
#' 
#' Test for qualitative interactions between treatment effects and patient subgroups.
#' Perform the testing based on the estimated treatment effects and their standard
#' errors. Output all the results related with qualitative interaction tests as 
#' a "qualint" object, which includes all the results related with the testing. 
#' Two common tests for qualitative interactions are included: IBGA and LRT, among
#' which IBGA is the default. Useful for three types of responses: continuous, binary
#' and survival data. Complete data is needed as input.
#' 
#' @param y response variable. A numeric vector for \code{type = "continuous"}. 
#' A numeric vector with values equal to either 1 or 0 for \code{type = "binary"}.
#' A "Surv" object for \code{type = "survival"} (The function \code{Surv()} in 
#' package \pkg{survival} produces such a matrix).
#' @param trtment treatment variable. A vector with different values representing
#' different treatment groups. Each element corresponds to the treatment the
#' patient received. Should have the same length as \code{y}.
#' @param subgrp patient subgroup variable. A vector with the same length as 
#' \code{y} and \code{trtment}, with different values representing different patient
#' subgroups. Each element corresponds to the patient subgroup the patient belongs
#' to.
#' @param type response type (see above). Three types of responses are included right
#' now: \code{continuous, binary, survival}. By default, it assumes
#' the response variable is continuous.
#' @param scale the scale type for treatment effects. When \code{type = "continuous"},
#' \code{scale = "MD"} by default (mean difference). When \code{type = "binary"}, three types
#' of scales are available, which are \code{scale = "RD"} (risk difference), 
#' \code{"RR"} (relative risk) or \code{"OR"} (odds ratio). 
#' When \code{type = "survival"}, \code{scale = "HR"} (hazard ratio). 
#' @param test testing method. Choose either \code{"IBGA"} (interval based graphical
#' approach) or \code{"LRT"} (Gail Simon likelihood ratio test).
#' @param alpha significance level. The type I error for qualitative 
#' interaction tesing. The default is 0.05.
#' @param plotout whether output the plot or not for \code{test = "IBGA"}.
#' There is no plot output for \code{test = "LRT"}.
#' 
#' @details
#' In order to test for qualitative interactions between treatment effects and patient
#' subgoups, estimated treatment effects and their standard errors are necessary. 
#' For continuous responses, mean difference is derived as the meansure of treatment
#' effects with with its standard error equal to \eqn{\sqrt{sd_1^2/n_1+sd_2^2/n_2}}.
#' For binary responses, three different scales are available to measure the treatment
#' effects: risk difference, log relative risk and log odds ratio. Their standard
#' errors could easily obtained according to formulas. For survival responses, the log
#' hazard ratio is used to evaluate the treatment effects. The cox regression model is
#' used in this function to estimate the log hazard ratio and also its standard error.
#' 
#' For the IBGA graph, however, we plot it according to common measures of treatment
#' effects instead of the one used in the calculation. For continuous responses, mean
#' difference is used since it is the common treament effect scale. For binary
#' responses, the function plots risk difference, relative risk, and odds ratio directly.
#' For survival responses, hazard ratios are plotted instead of log hazard ratios.
#' 
#' In the power calculation, this function assumes the estimated treatment effect scale
#' and its standard errors are equal to the true values. For IBGA method, an explicit
#' formula is available, so it is very easy to calculate the power. For LRT, a simulation
#' is used to assess the power since no explicit formula is available.   
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
#' \code{\link{plot.qualint}}, \code{\link{qualval}}
#' 
#' @examples
#' 
#' #### Continuous ####
#' ynorm <- rnorm(300)
#' trtment <- sample(c(0, 1), 300, prob = c(0.4, 0.6), 
#'                  replace = TRUE)
#' subgrp <- sample(c(0, 1, 2), 300, prob = c(1/3, 1/3, 1/3), 
#'                  replace = TRUE)
#' test1 <- qualint(ynorm, trtment, subgrp)
#' plot(test1)
#' print(test1)
#' coef(test1)
#' test2 <- qualint(ynorm, trtment, subgrp, plotout = TRUE)
#' test3 <- qualint(ynorm, trtment, subgrp, test = "LRT")
#' 
#' #### Binary ####
#' ybin <- sample(c(0, 1), 300, prob = c(0.3, 0.7), 
#'                replace = TRUE)
#' test4 <- qualint(ybin, trtment, subgrp, type = "binary")
#' test5 <- qualint(ybin, trtment, subgrp, type = "binary", 
#'                  scale = "RR", test = "LRT")
#' 
#' #### Survival ####
#' time <- rpois(300, 200)
#' censor <- sample(c(0, 1), 300, prob = c(0.7, 0.3), 
#'                  replace = TRUE)
#' test6 <- qualint(Surv(time, censor), trtment, subgrp)
#' test7 <- qualint(Surv(time, censor), trtment, subgrp, 
#'                  type = "survival", test = "LRT")
#' test8 <- qualint(Surv(time, censor), trtment, subgrp, 
#'                  test = "IBGA", plotout = TRUE)
#'                  
#' @importFrom survival Surv                
#' @export

#### Interval Based Graphical Approach ####

qualint <- function(y, trtment, subgrp, 
                    type = c("continuous", "binary", "survival"),
                    scale = c("MD", "RD", "RR", "OR", "HR"),
                    test = c("IBGA", "LRT"),
                    alpha = 0.05, plotout = FALSE) {
  
  type <- match.arg(type)
  if((class(y) == "Surv") & (type != "survival")) {
    type <- "survival"
    warning("y is a Surv object and type is set to survival")
  }
  
  test <- match.arg(test)
  
  scale <- match.arg(scale)
  if(type == "continuous") {  
    if(scale != "MD"){
      scale <- "MD"
      warning("scale is set to MD for continuous response")
    } 
  } else if(type == "survival") {
    if(scale != "HR") {
      scale <- "HR"
      warning("scale is set to HR for survival response")
    }
  } else if(type == "binary") {
    if(!any(c("RD", "RR", "OR") == scale)) {
      scale <- "RD"
      warning("scale is set to RD for binary response")
    } 
  }  
      
  alpha <- as.vector(alpha) 
  if(length(alpha) != 1) {
    warning("alpha has length > 1 and only the first element will be used")
    alpha <- alpha[1]
  }   
  if(is.null(alpha)) stop("argument alpha can not be NULL")
  alpha <- as.numeric(alpha)
  if(is.na(alpha) | (alpha >= 1) | (alpha <= 0))
    stop("alpha should be a numeric value between 0 and 1")
  
  this_call = match.call()
  
  trtment <- as.vector(trtment)
  subgrp <- as.vector(subgrp)  
  y <- drop(y)
  dimy <- dim(y)
  nrowy <- ifelse(is.null(dimy), length(y), dimy[1])
  lentrt <- length(trtment)
  lensbp <- length(subgrp)
  if((nrowy != lentrt) | (nrowy != lensbp) | (lentrt != lensbp))
    stop("numbers of observations in y, trtment, subgrp are not equal")
  
  subgrp <- ifelse(is.na(subgrp), NA, subgrp)
  subgrp <- factor(subgrp)
  sbp_ini <- levels(subgrp)
  sbp_lev <- sbp_ini
  nsbp <- length(sbp_lev)
    
  trtment <- ifelse(is.na(trtment), NA, trtment)
  trtment <- factor(trtment)
  trt_ini <- levels(trtment)
  trt_lev <- trt_ini
  ntrt <- length(trt_lev)
  if(ntrt != 2)
    stop("this function is only available for two treatments analysis")
  
  allout <- switch(type,
                   continuous = conout(y, trtment, trt_ini, subgrp, sbp_lev, nsbp),
                   binary = binout(y, trtment, trt_ini, subgrp, sbp_lev, nsbp, scale),
                   survival = surout(y, trtment, trt_ini, subgrp, sbp_lev, nsbp))
  betaout <- allout$beta
  
  testout <- switch(test, 
                    IBGA = ibgaout(betaout, nsbp, alpha),
                    LRT = lrtout(betaout, nsbp, alpha))
  
  CI_index <- qnorm(1 - alpha / 2)
  obj <- switch(test,
                IBGA = list(call = this_call, n = allout$n, type = type, alpha = alpha,  
                            treatment = allout$trtment, reference = allout$trtname, 
                            nsbp = nsbp, subgroup = sbp_lev, scale = scale,
                            effect = betaout[, 1], se = betaout[, 2],
                            LowerCI = betaout[, 1] - CI_index * betaout[, 2],
                            UpperCI = betaout[, 1] + CI_index * betaout[, 2],
                            test = test, index = testout$index, 
                            LowerTI = betaout[, 1] - testout$index * betaout[, 2],
                            UpperTI = betaout[, 1] + testout$index * betaout[, 2],
                            pvalue = testout$pvalue, power = testout$power, nobs = nrowy,
                            missing = which(is.na(subgrp) | is.na(subgrp) | is.na(y))),
                LRT = list(call = this_call, n = allout$n, type = type, alpha = alpha, 
                           treatment = allout$trtment, reference = allout$trtname, 
                           nsbp = nsbp, subgroup = sbp_lev, scale = scale,
                           effect = betaout[, 1], se = betaout[, 2], 
                           LowerCI = betaout[, 1] - CI_index * betaout[, 2],
                           UpperCI = betaout[, 1] + CI_index * betaout[, 2],
                           test = test, cvalue = testout$index, 
                           pvalue = testout$pvalue, power = testout$power, nobs = nrowy, 
                           missing = which(is.na(subgrp) | is.na(subgrp) | is.na(y))))
  
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




