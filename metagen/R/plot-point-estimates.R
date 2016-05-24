# Copyright (C) 2013-2014 Thomas W. D. Möbius (kontakt@thomasmoebius.de)
#
#     This program is free software: you can redistribute it and/or
#     modify it under the terms of the GNU General Public License as
#     published by the Free Software Foundation, either version 3 of the
#     License, or (at your option) any later version.
#
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
#     General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with this program. If not, see
#     <http://www.gnu.org/licenses/>.

################
### Boxplots ###
################

#' Plotting performance: Box plots for mean squared error
#'
#' Box plots for mean squared error.
#'
#' @param res The collected results from a computer experiment.
#' @return A plot object.
#' @export
boxMSE <- function(res) {
    boxp <- (  ggplot(res, aes(type, mse))
             + geom_boxplot(aes(fill=type))
             + geom_hline(yintercept=0, size=.42, linetype=4)
             + scale_x_discrete(expression(paste("Type of estimator")))
             + scale_y_continuous(name=expression("Mean squared error"))
             + coord_flip()
             + theme(  legend.position="none"
                     , axis.text = element_text(colour = "black"))
             )
   return(boxp)
}

#' Plotting performance: Box plots for bias
#'
#' Box plots for the bias.
#'
#' @param res The collected results from a computer experiment.
#' @return A plot object.
#' @export
boxBias <- function(res) {
    boxp <- (  ggplot(res, aes(type, bias))
             + geom_boxplot(aes(fill=type))
             + geom_hline(yintercept=0, size=.42, linetype=4)
             + coord_flip()
             + scale_x_discrete(expression(paste("Type of estimator")))
             + scale_y_continuous(name=expression("Bias"))
             + theme(  legend.position="none"
                     , axis.text = element_text(colour = "black"))
             )
   return(boxp)
}

#' Plotting performance: Box plots for standard deviation
#'
#' Box plots for standard deviation.
#'
#' @param res The collected results from a computer experiment.
#' @return A plot object.
#' @export
boxSD <- function(res) {
    boxp <- (  ggplot(res, aes(type, sqrt(variance)))
             + geom_boxplot(aes(fill=type))
             + geom_hline(yintercept=0, size=.3, linetype=4)
             + coord_flip()
             + scale_x_discrete(expression(paste("Type of estimator")))
             + scale_y_continuous(name=expression("Standard deviation"))
             + theme(  legend.position="none"
                     , axis.text = element_text(colour = "black"))
             )
   return(boxp)
}

##############################
### Scatterplots against h ###
##############################

#' Plotting performance: Scatter plots against heterogeneity
#'
#' Scatter plots of heterogeneity and mean squared error.
#'
#' @param res The collected results from a computer experiment.
#' @param ... further arguments to scale_y_continuous
#' @return A plot object.
#' @export
sctMSE <- function(res, ...) {
    p <- (  ggplot(res, aes(h, mse))
          + geom_point(aes(shape=type, colour=type))
          + geom_hline(yintercept=0, size=.3, linetype=4)
          + scale_x_continuous(expression(tau))
          + scale_y_continuous(expression("Mean squared error"), ...)
          + scale_shape_discrete(expression(paste("Type of ", tau, "-estimator")))
          + scale_colour_discrete(expression(paste("Type of ", tau, "-estimator")))
          + scale_linetype_discrete(expression("Estimated regression line"))
          + stat_smooth(  method="lm", se=FALSE, colour="black"
                        , aes(linetype=type))
          + theme(axis.text = element_text(colour = "black"))
          )
    return(p)
}

#' Plotting performance: Scatter plots against heterogeneity
#'
#' Scatter plots of heterogeneity and bias.
#'
#' @param res The collected results from a computer experiment.
#' @param ... further arguments to scale_y_continuous
#' @return A plot object.
#' @export
sctBias <- function(res, ...) {
    p <- (  ggplot(res, aes(h, bias))
          + geom_point(aes(shape=type, colour=type))
          + geom_hline(yintercept=0, size=.3, linetype=4)
          + scale_x_continuous(expression(tau))
          + scale_y_continuous(expression("Bias"), ...)
          + scale_shape_discrete(expression(paste("Type of ", tau, "-estimator")))
          + scale_colour_discrete(expression(paste("Type of ", tau, "-estimator")))
          + scale_linetype_discrete(expression("Estimated regression line"))
          + stat_smooth(method="lm", se=FALSE, colour="black",
                        aes(linetype=type))
          + theme(axis.text = element_text(colour = "black"))
          )
    return(p)
}

#' Plotting performance: Scatter plots against heterogeneity
#'
#' Scatter plots of heterogeneity and standard deviation.
#'
#' @param res The collected results from a computer experiment.
#' @param ... further arguments to scale_y_continuous
#' @return A plot object.
#' @export
sctSD <- function(res, ...) {
    p <- (  ggplot(res, aes(h, sqrt(variance)))
          + geom_point(aes(shape=type, colour=type))
          + geom_hline(yintercept=0, size=.3, linetype=4)
          + scale_x_continuous(expression(tau))
          + scale_y_continuous(expression("Standard deviation"), ...)
          + scale_shape_discrete(expression(paste("Type of ", tau,
                                                  "-estimator")))
          + scale_colour_discrete(expression(paste("Type of ", tau,
                                                   "-estimator")))
          + scale_linetype_discrete(expression("Estimated regression line"))
          + stat_smooth(method="lm", se=FALSE, colour="black",
                        aes(linetype=type))
          + theme(axis.text = element_text(colour = "black"))
          )
    return(p)
}
