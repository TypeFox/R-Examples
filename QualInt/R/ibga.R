#' Extract IBGA results from a "qualint" object
#' 
#' Extract the testing results of IBGA from a "qualint" object. Only available
#' when \code{test = "IBGA"}. Produce an error signal if \code{test = "LRT"}.
#' 
#' @param x a qualint object
#' 
#' @details 
#' This function extracts the results related with the IBGA graph from a "qualint"
#' object, which is a matrix which contains the estimating and testing information
#' of the treatment effects (see examples). FOr continuous responses, it is the 
#' mean difference. For binary responses, it is the risk difference, relative risk 
#' or odds ratio. For survival responses, it it hazard ratio.
#' 
#' @return A numeric matrix depending on \code{scale} (see above).
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
#' @seealso \code{\link{coef.qualint}}
#' 
#' @examples
#' ynorm <- rnorm(300)
#' trtment <- sample(c(0, 1), 300, prob = c(0.4, 0.6), 
#'                   replace = TRUE)
#' subgrp <- sample(c(0, 1, 2), 300, prob = c(1/3, 1/3, 1/3), 
#'                  replace = TRUE)
#' test1 <- qualint(ynorm, trtment, subgrp)
#' ibga(test1)
#' 
#' @export

ibga <- function(x){
  
  if(class(x) != "qualint")
    stop("this function is only available for qualint object")
  
  if(x$test == "LRT")
    stop("this function is not available for LRT")
  
  if((x$type == "continuous") | ((x$type == "binary") & (x$scale == "RD")))
    out <- matrix(c(x$effect, x$se, x$LowerCI, x$UpperCI, x$LowerTI, x$UpperTI), 
                  x$nsbp, 6, dimnames = list(x$subgroup, c("Estimate", "Std. Error", 
                                                           "Lower CI", "Upper CI",
                                                           "Lower TI", "Upper TI")))
  else out <- matrix(c(exp(x$effect), abs(exp(x$effect) * x$se), exp(x$LowerCI), 
                       exp(x$UpperCI), exp(x$LowerTI), exp(x$UpperTI)), 
                     x$nsbp, 6, dimnames = list(x$subgroup, c("Estimate", "Std. Error", 
                                                              "Lower CI", "Upper CI",
                                                              "Lower TI", "Upper TI")))
  
  class(out) <- c("ibga", as.character(x$scale), class(out))
  out
  
}

#' @export
print.ibga <- function(x, digits = max(3, getOption("digits") - 3), ...){
  switch(class(x)[2],
         MD = cat("IBGA Results for Mean Difference:\n"),
         RD = cat("IBGA Results for Risk Difference:\n"),
         RR = cat("IBGA Results for Relative Risk:\n"),
         OR = cat("IBGA Results for Odds Ratio:\n"),
         HR = cat("IBGA Results for Hazard Ratio:\n"))
  class(x) <- "matrix"
  print(x, digits = digits)
}




