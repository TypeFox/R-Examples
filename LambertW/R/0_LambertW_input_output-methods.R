#' @title Methods for Lambert W input and output objects
#' @name LambertW_input_output-methods
#' @aliases print.LambertW_input plot.LambertW_input print.LambertW_output plot.LambertW_output
#' 
#' @description
#' S3 methods for Lambert W input and output objects 
#' (created by \code{\link{create_LambertW_input}} and \code{\link{create_LambertW_output}}).
#' 
#' @param x object of class \code{LambertW_input} or \code{LambertW_output}.
#' @param xlim lower and upper limit of x-axis for cdf and pdf plots. If \code{NULL},
#' it tries to determine good limits based on the family type of the distribution
#' and the quantiles.  Most of the times it will show the pdf and cdf from the 0.5\% to
#' 99.5\% quantile.
#' @param \dots further arguments passed to or from other methods.
#' @keywords hplot print
#' @examples
#' 
#' # create a Normal(1, 2) input
#' Gauss.input <- create_LambertW_input("normal", beta = c(1, 2))
#' plot(Gauss.input)
#' # make it a bit heavy tailed (beta in theta comes from Gauss.input)
#' LW.Gauss <- create_LambertW_output(LambertW.input = Gauss.input, 
#'                                    theta = list(delta = c(0.3)))
#' LW.Gauss  # print a nice overview in the console
#' plot(LW.Gauss)
#' 
#' # draw random sample
#' LW.Gauss$r(n=10)
#' Gauss.input$r(n=10)
#' # quantiles
#' LW.Gauss$q(p=0.6)
#' Gauss.input$q(p=0.6)
#' 
NULL