#' @name CALF-package
#' @aliases CALF-package
#' @title Coarse Approximation Linear Function
#' @description Forward selection linear regression greedy algorithm.
#' @author {Stephanie Lane [aut, cre],\cr
#'    Clark Jeffries [aut], \cr
#'    Diana Perkins [aut], cr
#' }
#' Maintainer: Stephanie Lane \email{slane@@unc.edu}
#' @importFrom stats t.test
#' @import ggplot2
#' @keywords calf
#' @details The Coarse Approximation Linear Function (CALF) algorithm is a type of forward selection
#' linear regression greedy algorithm. Nonzero weights are restricted to the values +1 and -1.
#' The number of nonzero weights used is limited by a parameter. Samples are controls (at least 2) and cases (at least 2).
#' A data matrix consists of a distinguished column that labels every row as either a control (0) or a case (1).
#' Other columns (at least one) contain real number marker measurement data.
#' Another input is a limit (positive integer) on the number of markers that can be selected for use in a linear sum.
#' The present version uses as a score of differentiation the two-tailed, two sample unequal variance Student t-test p-value.
#' Thus, any real-valued function applied to all samples generates values for controls and cases that are used to calculate the score.
#' CALF selects the one marker (first in case of tie) that best distinguishes controls from cases (score is smallest p-value).
#' CALF then checks the limit. If the number of selected markers is the limit, CALF ends.
#' Else, CALF seeks a second marker, if any, that best improves the score of the sum function generated
#' by adding the newly selected marker to the previous markers with weight +1 or weight -1.
#' The process continues until the limit is reached or until no additional marker can be included in the sum to improve the score.
NULL
