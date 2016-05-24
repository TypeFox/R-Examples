#' OrthoPanels: Orthogonalized Panel Model
#'
#' This package includes the function \code{opm()}, which implements the
#' orthogonal reparameterization approach recommended by Lancaster
#' (2002) to estimate dynamic panel models with fixed effects (and
#' optionally: panel specific intercepts). The OLS estimator for such
#' models is biased with a fixed (small) \eqn{N} (Nickell 1981). Equivalently,
#' a maximum likelihood estimation leads to an incidental parameters
#' problem (Neyman and Scott 1948; Lancaster 2000). The approach by
#' Lancaster (2002) uses an orthogonal reparameterization of the fixed
#' effects to produce a likelihood-based estimator of the remaining
#' parameters that is exact and consistent as \eqn{N} approaches infinity
#' for \eqn{T} greater than or equal to 2.
#'
#' @references
#'
#' Lancaster, T. (2000) The incidental parameter problem since 1948.
#' \emph{Journal of Econometrics}, \bold{95}, 391--413.
#' 
#' Lancaster, T. (2002) Orthogonal parameters and panel data.
#' \emph{Review of Economic Studies}, \bold{69}, 647--666.
#'
#' Neyman, J. and Scott, E. L. (1948) Consistent estimation from
#' partially consistent observations. \emph{Econometrica}, \bold{16},
#' 1--32.
#' 
#' Nickell, S. (1981) Biases in dynamic models with fixed effects.
#' \emph{Econometrica}, \bold{49}, 1417--1426.
#'
#' @name OrthoPanels-package
#' @docType package
NULL

#' Responses from the 2010 British Election Study
#'
#' A survey of 1845 respondents using 3 waves of panel survey data
#' from the 2010 British Election Study. The variables are as follows:
#'
#' \itemize{
#'   \item n case number
#'   \item t time wave
#'   \item Econ Assessment of change in the national economic situation over the past 12 months
#'         (1-5, 1=\sQuote{got a lot worse}, 5=\sQuote{got a lot better})
#'   \item Clegg Evaluation of Liberal Party leader Nick Clegg
#'         (0-10, 0=\sQuote{strongly dislike} and 10=\sQuote{strongly like})
#'   \item Brown Evaluation of Labour Party leader Gordon Brown
#'   \item Cameron Evaluation of Conservative Party leader David Cameron
#'   \item Approve Approval of the government, as expressed by feeling about the ruling Labour Party
#'         (0-10, 0=\sQuote{strongly dislike}, 10=\sQuote{strongly like})
#'   \item NHS Assesment of the current government's handling of the National Health
#'         Service (1-5, 1=\sQuote{very badly}, 5=\sQuote{very well})
#'   \item Terror Assesment of the current government's handling of terrorism
#'         (1-5, 1=\sQuote{very badly}, 5=\sQuote{very well})
#'   \item PID Personal identification with the Labour Party
#'         (0/1, 0=\sQuote{no}, 1=\sQuote{yes})
#'   \item Tax Preference for policy on taxes and health and social spending
#'         (0-10, 0=\sQuote{cut taxes a lot and spend much less},
#'         10=\sQuote{increase taxes a lot and spend much more})
#' }
#'
#' @format A data frame with 5535 rows and 11 variables
#' @name BES_panel
NULL
