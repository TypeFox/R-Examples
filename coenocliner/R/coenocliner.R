##' @title A coenocline simulation package for R
##'
##' @description \pkg{coenocliner} provides a simple, easy interface for simulating species abundance (counts) or occurrence along gradients.
##'
##' @details One of the key ways quantitative ecologists attempt to understand the properties and behaviour of the methods they use or dream up is through the use of simulated data. \pkg{coenocliner} is an R package that provides a simple interface to coenocline simulation.
##'
##' Species data can be simulated from a number of species response models
##' \itemize{
##'   \item Gaussian response
##'   \item Generalised Beta response
##' }
##' and random count or occurrence data can be simulated from suitably parameterised response models by using the output from the response model as the mean or expectation of one of a number of statistical distributions
##' \itemize{
##'   \item Poisson
##'   \item Negative Binomial
##'   \item Bernoulli
##'   \item Binomial
##'   \item Beta-Binomial
##'   \item Zero-inflated Poisson (ZIP)
##'   \item Zero-inflated Negative Binomial (ZINB)
##'   \item Zero-inflated Binomial (ZIB)
##'   \item Zero-inflated Beta-Binomial (ZIBB)
##' }
##' from which random draws are made.
##'
##' @author Gavin L. Simpson
##' @docType package
##' @name coenocliner
##' @aliases coenocliner-package
##' @keywords package
##' @seealso \code{\link{coenocline}} for simulating species data, \code{\link{distributions}} for details of the error distributions tha can be used for simulations, and \code{\link{species-response}} for details on the available species response models and the parameters required to use them.
NULL
