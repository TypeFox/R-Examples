#' Extract estimation results from a "qualint" object
#' 
#' Similar to other coef methods, this function 
#' extracts the estimation results from an "qualint" object.
#' 
#' @param object a qualint object
#' @param ... not used. Additional print arguments
#' 
#' @details
#' This function extracts the results related with estimating results of the
#' interaction test from a "qualint" object. It returns a matrix which contains 
#' the estimating results of the treatment effects which will used directly in
#' the testing process later. FOr continuous responses, it is the 
#' mean difference. For binary responses, it is the risk difference, log relative 
#' risk or log odds ratio. For survival responses, it it the log hazard ratio.
#' 
#' @return A numeric matrix (see above). For \code{type = "continuous"},
#' it is the estimation results for mean difference. For \code{type = "binary"},
#' it is for risk difference, log relative risk or log odds ration. For
#' \code{type = "survival"}, it is for log hazard ratio.
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
#' @seealso \code{\link{print.qualint}}, \code{\link{plot.qualint}}
#' 
#' @examples
#' ynorm <- rnorm(300)
#' trtment <- sample(c(0, 1), 300, prob = c(0.4, 0.6), 
#'                   replace = TRUE)
#' subgrp <- sample(c(0, 1, 2), 300, prob = c(1/3, 1/3, 1/3), 
#'                  replace = TRUE)
#' test1 <- qualint(ynorm, trtment, subgrp)
#' coef(test1)
#' 
#' @export

coef.qualint <- function(object, ...){
  
  outcoef <- matrix(c(object$effect, object$se, object$LowerCI, 
                      object$UpperCI), object$nsbp, 4,
                    dimnames = list(object$subgroup, 
                                    c("Estimate", "Std. Error", 
                                      "Lower CI", "Upper CI")))
  class(outcoef) <- c("coef", as.character(object$scale), class(outcoef))
  outcoef
  
}

#' @export
print.coef <- function(x, digits = max(3, getOption("digits") - 3), ...){
  switch(class(x)[2],
         MD = cat("Estimation Results for Mean Difference:\n"),
         RD = cat("Estimation Results for Risk Difference:\n"),
         RR = cat("Estimation Results for log(Relative Risk):\n"),
         OR = cat("Estimation Results for log(Odds Ratio):\n"),
         HR = cat("Estimation Results for log(Hazard Ratio):\n"))
  class(x) <- "matrix"
  print(x, digits = digits)
}
