#' intercure: Cure rate regression models for interval censored data.
#'
#' The intercure package provides two main functions: \code{inter_bch} and
#' \code{inter_frailty}. These are essentialy algorithms for estimating the cure
#' fraction with promotion time and frailty model, respectively. The
#' \code{inter_frailty_cl} function provides support for analysing clustered
#' datasets using the frailty model. For generating datasets based on these two
#' models, the package provides the \code{sim_bch}, \code{sim_frailty}
#' and \code{sim_frailty_cl} functions, the last providing clustered datasets.
#'
#'
#' @section intercure functions:
#' \code{inter_bch}, \code{inter_frailty}, \code{inter_frailty_cl}, \code{sim_bch}, \code{sim_frailty}, \code{sim_frailty_cl}
#'
#' @docType package
#' @name intercure
#' @import foreach
NULL
