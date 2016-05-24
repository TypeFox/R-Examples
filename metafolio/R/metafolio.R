#' metafolio: An R package to simulate metapopulations for
#' portfolio optimization
#'
#' The \pkg{metafolio} \R package is a tool to simulate metapopulations and
#' apply financial portfolio optimization concepts. The package was originally
#' written for salmon simulations, so some of the language refers to
#' salmon-specific terminology, but the package could be used and/or adopted
#' for other taxonomic groups.
#'
#' The main simulation function is \code{\link{meta_sim}}. This function takes
#' care of running an individual simulation iteration. The package also
#' contains functions for exploring conservation scenarios with these
#' simulations (see the "Assessing multiple conservation scenarios" section
#' below), and find optimal conservation strategies (see the "Portfolio
#' optimization section" below).
#'
#' @section Running a simulation once:
#'
#' To run a single simulation iteration, see the function
#' \code{\link{meta_sim}}. To plot the output from one of these simulations, see
#' the function \code{\link{plot_sim_ts}}.
#'
#' @section Assessing multiple conservation scenarios:
#'
#' You can use \code{\link{run_cons_plans}} to run \code{\link{meta_sim}} for
#' multiple iterations and across multiple conservation strategies. These
#' strategies could focus on the spatial distribution of conservation or on the
#' number of populations conserved.
#'
#' The function \code{\link{plot_cons_plans}} can plot the output from
#' \code{\link{run_cons_plans}}.
#'
#' @section Specifying environmental patterns:
#'
#' When you run \code{\link{meta_sim}} you can specify the environmental signal.
#' One of the arguments is a list of options to pass to
#' \code{\link{generate_env_ts}}, which controls the environmental pattern.
#'
#' @section Diagnostic plots:
#'
#' \pkg{metafolio} contains some additional plotting functions to inspect the
#' spawner-return relationships and the correlation between returns:
#' \code{\link{plot_rickers}}, and
#' \code{\link{plot_correlation_between_returns}}.
#'
#' @section Portfolio optimization:
#'
#' \pkg{metafolio} also contains some experimental functions for finding optimal
#' conservation strategies (an efficient frontier). This is analogous to
#' financial portfolio where the goal is to find the investment weights that
#' maximizes expected return for a level of expected risk, or vice-versa.
#' Presently, these functions rely on Monte Carlo sampling, and so are rather
#' slow.
#'
#' For this purpose, the function \code{\link{create_asset_weights}} can
#' generate a matrix of asset weights, which can then be passed to
#' \code{\link{monte_carlo_portfolios}} to do the optimization itself.
#' \code{\link{plot_efficient_portfolios}} can be used to plot the optimization
#' output.
#'
#' See the package vignette \code{vignette("metafolio")} for more extensive
#' explanation of how to use \pkg{metafolio} along with some examples.
#'
#' @docType package
#' @name metafolio
#' @importFrom Rcpp cppFunction
NULL
