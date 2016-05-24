

#' Plotting longitudinal series using nicer graphical parameters.
#' 
#' The aim is to create nicer longitudinal series plots.
#' 
#' \tabular{ll}{ Package: \tab Eplot\cr Type: \tab Package\cr Version: \tab
#' 1.0\cr Date: \tab 2014-07-30\cr License: \tab GPL-2 \cr } The idea is to
#' facilitate the use of nicer graphical parameters. The user have the choice
#' to keep the new set of graphical parameters or to revert to her initial one.
#' Other functions include multivariate plot, plot with vertical axis on both
#' left and right hand sides, and plot which superimpose prediction intervals
#' from an AR-ARCH model.
#' 
#' @name Eplot
#' @aliases Eplot-package Eplot
#' @docType package
#' @author Eran Raviv
#' 
#' Maintainer: Eran Raviv \email{Eran.Raviv@apg-am.nl}
#' @examples
#' 
#'  par(mfrow = c(2,1))
#'  out <- FCIplot(rnorm(100),plott=T,k=40)
#'  plott(out,tit="The out-of-sample standard deviation") 




