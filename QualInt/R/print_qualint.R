#' Print a summary of a "qualint" object
#' 
#' Similar to other print methods, this function prints a summary from an "qualint"
#' object.
#' 
#' @param x a qualint object
#' @param digits significant digits in printout
#' @param ... not used. Additional print arguments
#' 
#' @details
#' This function prints the important information in a format that's easier for
#' people to understand from a "qualint" object (see examples). 
#' 
#' @return A summary of the testing results (see above).
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
#' @seealso \code{\link{coef.qualint}}, \code{\link{plot.qualint}}
#' 
#' @examples
#' 
#' ynorm <- rnorm(300)
#' trtment <- sample(c(0, 1), 300, prob = c(0.4, 0.6), 
#'                   replace = TRUE)
#' subgrp <- sample(c(0, 1, 2), 300, prob = c(1/3, 1/3, 1/3), 
#'                  replace = TRUE)
#' test1 <- qualint(ynorm, trtment, subgrp)
#' print(test1)
#' 
#' @export

print.qualint <- function(x, digits = max(3, getOption("digits") - 3), ...){
  
  if(is.null(x$reference))
    dimnam <- c("reference", "treatment")
  else dimnam <- c(paste("reference", x$reference, sep = ":"),
                   paste("treatment", 
                         x$treatment[x$treatment != x$reference], 
                         sep = ":"))
  
  cat("\nCall:\n")
  print(x$call)
  
  cat("\nType:\n")
  cat(x$type)
  
  #if(x$type == "continuous") {
  #  cat("\n\nMean:\n")
  #  print(matrix(c(x$mean[[1]], x$mean[[2]], x$effect),
  #               nrow = 3, byrow = TRUE,
  #               dimnames = list(c(dimnam, "difference"), x$subgroup)), 
  #        digits = digits)
  #} else if(x$type == "binary") {
  #  cat("\n\nEvent Risk:\n")
  #  print(matrix(c(x$risk[[1]], x$risk[[2]], switch(x$scale,
  #                                                  RD = x$effect,
  #                                                  RR = exp(x$effect),
  #                                                  OR = exp(x$effect))),
  #               nrow = 3, byrow = TRUE,
  #               dimnames = list(c(dimnam, 
  #                                 switch(x$scale, 
  #                                        RD = "risk diff", 
  #                                        RR = "relative risk",
  #                                        OR = "odds ratio")), 
  #                               x$subgroup)), 
  #        digits = digits)
  #} else if(x$type == "survival") {
  #  cat("\n\nEvent Rate:\n")
  #  print(matrix(c(x$rate[[1]], x$rate[[2]], exp(x$effect)),
  #               nrow = 3, byrow = TRUE,
  #               dimnames = list(c(dimnam, "hazard ratio"), x$subgroup)), 
  #        digits = digits)
  #}
  
  #cat("\n\nSample Size:\n")
  #print(matrix(c(x$n[[1]], x$n[[2]]),
  #             nrow = 2, byrow = TRUE,
  #             dimnames = list(dimnam, x$subgroup)))
  
  if(x$type == "continuous") cat("\n\nEstimating Results for Mean Difference:\n")
  else if(x$type == "binary") switch(x$scale, 
                                     RD = cat("\n\nEstimating Results for Risk Difference:\n"),
                                     RR = cat("\n\nEstimating Results for log(Relative Risk):\n"),
                                     OR = cat("\n\nEstimating Results for log(Hazard Ratio):\n"))
  else if(x$type == "survival") cat("\n\nEstimating Results for log(Hazard Ratio):\n")
  else if(x$type == "unknown") cat("\n\nEstimating Results for treatment effects:\n")
  print(matrix(c(x$effect, x$se, x$LowerCI, x$UpperCI), 
               x$nsbp, 4,
               dimnames = list(x$subgroup, c("Estimate", "Std. Error", 
                                             "Lower CI", "Upper CI"))),
        digits = digits)
  
  cat("\nTest:\n")
  cat(x$test)
  
  cat("\n\np-value:\n")
  cat(format(x$pvalue, digits = digits))
  
  cat("\n\nPower:\n")
  cat(format(x$power, digits = digits))
  
  cat("\n\nAlpha:\n")
  cat(x$alpha)
  cat("\n\n")
  
}
