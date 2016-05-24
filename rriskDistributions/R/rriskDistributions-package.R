
################################################################################
################################################################################
#' This packages provides a collection of functions for estimation parameters of
#' continuous or discrete distributions (related to the \code{rrisk} project)
#' to given data or to known quantiles.
#'
#' This package is a part of the \code{rrisk} project and contains functions for
#' fitting distributions to given data or by known quantiles. This package does
#' not depend on the whole \code{rrisk} project and can be used separately. The
#' \code{rrisk} project can be downloaded from \url{http://www.bfr.bund.de/cd/52158}.
#' \cr \cr
#' The main functions \code{fit.perc} and \code{fit.cont} call a GUI that allows users
#' to choose an appropriate distribution family to given data or to known quantiles
#' without any knowledge of the \acronym{R} syntax.
#'
#' @name rriskDistributions-package
#' @aliases rriskDistributions
#' @docType package
#' @concept Beta, Normal, Lognormal, Cauhy, Chi-Quadrat, Logistic, Student's t,
#' Exponential, F, Gamma, Weibull, rrisk, stat-up
#' @title Fitting distributions to given data or known quantiles
#' @note Fitting by given quantiles: a typical application is the definition of a distribution based on expert
#' opinion on some quantiles (e.g., the 2.5th, median and 97.5th) of the trial
#' to be modelled. \code{rrisk} has a functionality, to fit all continuous or
#' discrete distributions simultaneously without urging the user to specify the
#' distribution family in advance.
#' @keywords package
#' @author Natalia Belgorodski \email{belgorodski@@stat-up.de} (\acronym{STAT-UP} Statistical Consulting), \cr
#' Matthias Greiner \email{matthias.greiner@@bfr.bund.de} (Federal Institute for Risk Assessment, Germany), \cr
#' Kristin Tolksdorf \email{kristin.tolksdorf@@bfr.bund.de} (Federal Institute for Risk Assessment, Germany), \cr
#' Katharina Schueller \email{schueller@@stat-up.de} (\acronym{STAT-UP} Statistical Consulting)
#' @examples
#' q <- stats::qweibull(p = c(0.025, 0.5, 0.975), shape = 2, scale = 3)
#' get.weibull.par(q = q)
#' q <- stats::qweibull(p = c(0.025, 0.5, 0.975), shape = 0.01, scale = 1)
#' get.weibull.par(q = q)
#' 
#' p <- c(0.025, 0.50, 0.975)
#' q <- c(9.68, 29.2, 50.98)
#' fit.results <- rriskFitdist.perc(p, q, show.output = FALSE)
#' plotDiagnostics.perc(fit.results)
#'
#' p <- c(0.25, 0.50, 0.75)
#' q <- c(9.68, 29.2, 50.98)
#' fit.results <- rriskFitdist.perc(p, q, show.output = FALSE)
#' plotDiagnostics.perc(fit.results)
#' plotDiagnostics.perc(fit.results, tolPlot = 2)
#' 
#' \dontrun{
#'   if( class(tcltk::tclRequire("Tktable")) == "tclObj" ) {
#'     res.fitcont <- fit.cont(data2fit = rnorm(100))
#'     res.fitcont
#'   }
#'   if( class(tcltk::tclRequire("Tktable")) == "tclObj" ) {
#'     res.fitperc <- fit.perc()
#'     res.fitperc
#'   }
#' }
#' 
NULL
