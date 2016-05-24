#' Mann-Whitney-Wilcoxon ranks test when data are in groups.
#'
#' Calculate a Mann-Whitney-Wilcoxon test for a difference between treatment
#' levels using nested ranks.  This test can be used when observations are
#' structured into several groups and each group has received both treatment
#' levels.  The p-value is determined via bootstrapping.  This test is intended
#' to be analogous to a mixed-model extension of the \code{\link{wilcox.test}},
#' for which treatment is a fixed effect and group membership is a random
#' effect.
#'
#' The main function is \code{\link{nestedRanksTest}}, which includes a
#' formula interface implementing the familiar \code{"|"} syntax for
#' specifying group membership on the right-hand side of the formula.
#' The value returned is a list of class \code{'htest_boot'}, which
#' extends class \code{'htest'}.  \code{print} and \code{plot} methods
#' are provided to print and visualise results.
#'
#' These statistical tools were developed in collaboration with Peter E.
#' Smouse (Rutgers University) and Victoria L. Sork (UCLA) and were funded
#' in part by U.S. National Science Foundation awards NSF-DEB-0514956 and
#' NSF-DEB-0516529.
#'
#' @references
#' Thompson, P. G., Smouse, P. E., Scofield, D. G. and Sork, V. L. (2014)
#' What seeds tell us about birds: a multi-year analysis of acorn woodpecker
#' foraging movements.  \emph{Movement Ecology} 2:12.
#' \url{http://www.movementecologyjournal.com/content/2/1/12}
#'
#' \url{https://github.com/douglasgscofield/nestedRanksTest}
#'
#' @docType package
#'
#' @name nestedRanksTest-package
#'
NULL
