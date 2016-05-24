#' R-package for qualitative interaction test
#' 
#' Test for qualitative interactions between treatment effects and 
#' patient subgroups for continuous, binary and suvival responses. The term treatment 
#' effect means a comparison result between two treatments within each patient 
#' subgroup.
#' 
#' @details
#' 
#' \tabular{ll}{
#' Package: \tab QualInt \cr
#' Type:    \tab Package \cr
#' Version: \tab 1.0.0 \cr
#' Date:    \tab 2014-10-13 \cr
#' License: \tab GPL-2 
#' }
#' 
#' This package could be used to calculate the pvalue and power for qualitative
#' interaction testing. Two testing methods are included in the package, which are
#' Interval Based Graphical Approach and Gail Simon Likelihood Ratio Test.
#' 
#' For a complete list of all the functions available in this package with
#' individual help pages, use library(help = "QualInt")
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
#' @examples
#' 
#' test9 <- qualval(effect = c(1.0, 0.5, -2.0), se = c(0.86, 0.64, 0.32))
#' print(test9)
#' plot(test9)
#' 
#' #### Continuous ####
#' ynorm <- rnorm(300)
#' trtment <- sample(c(0, 1), 300, prob = c(0.4, 0.6), 
#'                   replace = TRUE)
#' subgrp <- sample(c(0, 1, 2), 300, prob = c(1/3, 1/3, 1/3), 
#'                  replace = TRUE)
#' test1 <- qualint(ynorm, trtment, subgrp)
#' test2 <- qualint(ynorm, trtment, subgrp, test = "LRT")
#' plot(test1)
#' print(test1)
#' coef(test1)
#' ibga(test1)
#' 
#' #### Binary ####
#' ybin <- sample(c(0, 1), 300, prob = c(0.3, 0.7), 
#'                replace = TRUE)
#' test4 <- qualint(ybin, trtment, subgrp, type = "binary")
#' 
#' #### Survival ####
#' time <- rpois(300, 200)
#' censor <- sample(c(0, 1), 300, prob = c(0.7, 0.3), 
#'                  replace = TRUE)
#' test6 <- qualint(Surv(time, censor), trtment, subgrp)
#' 
#' @docType package
#' 
#' @name QualInt-package

NULL