# Copyright (C) 2014 Thomas W. D. Möbius (kontakt@thomasmoebius.de)
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

turn2perc <- function(confidence) {
     return(factor(  paste(confidence*100, "%")
                   , ordered=TRUE))
}

#' Plot pivots: Interval estimates of the heterogeneity
#'
#' @param cnfh interval estimates of the heterogeneity.
#' @export
plotHeterogeneityInterval <- function(cnfh) {
    cnfh$confidence <- turn2perc(cnfh$confidence)
    return(ggplot(cnfh, aes(type, 0, ymin = lower, ymax = upper))
           + geom_linerange(size = 3.42)
           + scale_y_continuous(expression(tau))
           + scale_x_discrete(name = "")
           + coord_flip()
           + facet_wrap(~confidence, ncol=1)
           + theme(axis.text = element_text(colour = "black"))
           )
}

#' Plot pivots: Interval estimates of the heterogeneity
#'
#' @param cnfr interval estimates of the heterogeneity.
#' @export
plotCoefficientInterval <- function(cnfr) {
    cnfr$confidence <- turn2perc(cnfr$confidence)
    return( ggplot(cnfr, aes(method, 0, ymin = lower, ymax = upper))
           + geom_linerange(size = 3.42)
           + scale_y_continuous(expression(beta[2]))
           + scale_x_discrete(name = "")
           + geom_hline(yintercept = 0, size = 0.3, linetype=4)
           + coord_flip()
           + facet_wrap(~confidence, ncol=1)
           + theme(axis.text = element_text(colour = "black"))
           )
}

################
### Boxplots ###
################

#' Plotting performance: Box plots for target value confidence-coverage
#'
#' @param res The collected results from a computer experiment.
#' @return A plot object.
#' @export
boxByConfidence <- function(res) {
    boxp <- (  ggplot(res, aes(factor(confidence, ordered=T), confidence-coverage))
             + geom_boxplot(aes(fill=type))
             + geom_hline(yintercept=0, size=.42, linetype=4)
             + coord_flip()
             + scale_fill_discrete(expression("Type of method\nused for\ninterval estimation"))
             + scale_x_discrete(expression(paste("Aspired confidence level")))
             + ylab(expression("Confidence level - Coverage probability"))
             + theme(  legend.position="bottom"
                     , axis.text = element_text(colour = "black"))
             )
   return(boxp)
}

#' Plotting performance: Box plots for target value confidence-coverage
#'
#' @param res The collected results from a computer experiment.
#' @return A plot object.
#' @export
boxByType <- function(res) {
    boxp <- (  ggplot(res, aes(type, confidence-coverage))
             + geom_boxplot(aes(fill=factor(confidence, ordered=T)))
             + geom_hline(yintercept=0, size=.42, linetype=4)
             + coord_flip()
             + scale_x_discrete(expression(paste("Method used for interval estimation")))
             + scale_fill_discrete(expression(paste("Aspired\nconfidence\nlevel")))
             + ylab(expression("Confidence level - Coverage probability"))
             + theme(  legend.position="bottom"
                     , axis.text = element_text(colour = "black"))
             )
   return(boxp)
}

#' Plotting performance: Box plots for target value confidence-coverage
#'
#' @param res The collected results from a computer experiment.
#' @return A plot object.
#' @export
boxByMethod <- function(res) {
    boxp <- (  ggplot(res, aes(type, confidence-coverage))
             + geom_boxplot(aes(fill=method))
             + geom_hline(yintercept=0, size=.42, linetype=4)
             + coord_flip()
             + scale_x_discrete(expression("Underlying type used for interval estimation"))
             + scale_fill_discrete(expression("Type of interval estimate"))
             + ylab(expression("Confidence level - Coverage probability"))
             + theme(  legend.position="bottom"
                     , axis.text = element_text(colour = "black"))
             )
   return(boxp)
}

##############################################################
### Scatter plots against h and propagated confidence level ###
##############################################################

#' Plotting performance: Scatter plot against heterogeneity
#'
#' @param res The collected interval results from a computer experiment.
#' @return A plot object.
#' @export
sctVersusH <- function(res) {
    res$colour <- turn2perc(res$confidence)
    p <- (  ggplot(res, aes(h, confidence-coverage))
          + geom_point(aes(colour=colour), alpha=.729)
          + scale_colour_discrete(expression("Aspired confidence level"))
          + geom_hline(yintercept=0, size=.42, linetype=3)
          + scale_x_continuous(expression(tau))
          + ylab(expression("Confidence level - Coverage probability"))
          + stat_smooth(method="lm", colour="black")
          + theme(  legend.position="bottom"
                  , axis.text = element_text(colour = "black"))
          )
    return(p)
}

#' Plotting performance: Scatter plot against heterogeneity
#'
#' @param res The collected results from a computer experiment.
#' @return A plot object.
#' @export
sctVersusC <- function(res) {
    p <- (  ggplot(res, aes(confidence, confidence-coverage))
          + geom_point(aes(colour=h), alpha=.729)
          + scale_colour_gradient(expression("True heterogeneity")
                                  , low=muted("blue"), high=muted("red")
                                  , space="Lab"
                                  )
          + geom_hline(yintercept=0, size=.42, linetype=3)
          + scale_x_continuous(expression(tau))
          + ylab(expression("Confidence level - Coverage probability"))
          + stat_smooth(method="lm", colour="black") # aes(colour=type))
          + theme(  legend.position="bottom"
                  , axis.text = element_text(colour = "black"))
          )
    return(p)
}

###################################
### Scatterplots against d_mean ###
###################################

#' Plotting performance: Scatter plot against heterogeneity
#'
#' @param res The collected results from a computer experiment.
#' @return A plot object.
#' @export
sdmByType <- function(res) {
    sctp <- (  ggplot(res, aes(d_mean, confidence-coverage))
             + geom_point(alpha=5/8, aes(colour=type))
             + geom_hline(yintercept=0, size=.3)
             + scale_x_continuous(expression(bar(delta)))
             + ylab(expression("Confidence level - Coverage probability"))
             + estColourPalette
             + facet_grid(confidence ~ method))
    return(sctp)
}

#' Plotting performance: Scatter plot against heterogeneity
#'
#' @param res The collected results from a computer experiment.
#' @return A plot object.
#' @export
sdmByMethod <- function(res) {
    sdmp <- (  ggplot(res, aes(d_mean, confidence-coverage))
             + geom_point(alpha=5/8, aes(colour=factor(confidence)))
             + geom_hline(yintercept=0, size=.3)
             + scale_x_continuous(expression(bar(delta)))
             + ylab(expression("Confidence level - Coverage probability"))
             + facet_grid(. ~ method)
             + cnfColourPalette)
    return(sdmp)
}

#################################
### Scatterplots against d_sd ###
#################################

#' Plotting performance: Scatter plot against heteroscedasticity
#'
#' @param res The collected results from a computer experiment.
#' @return A plot object.
#' @export
sdsByType <- function(res) {
    sctp <- (  ggplot(res, aes(d_sd, confidence-coverage))
             + geom_point(alpha=5/8, aes(colour=type))
             + geom_hline(yintercept=0, size=.3)
             + scale_x_continuous(expression(sd(delta)))
             + ylab(expression("Confidence level - Coverage probability"))
             + estColourPalette
             + facet_grid(confidence ~ method))
    return(sctp)
}

#' Plotting performance: Scatter plot against heteroscedasticity
#'
#' @param res The collected results from a computer experiment.
#' @return A plot object.
#' @export
sdsByMethod <- function(res) {
    sdsp <- (  ggplot(res, aes(d_sd, confidence-coverage))
             + geom_point(alpha=5/8, aes(colour=factor(confidence)))
             + geom_hline(yintercept=0, size=.3)
             + scale_x_continuous(expression(sd(delta)))
             + ylab(expression("Confidence level - Coverage probability"))
             + facet_grid(. ~ method)
             + cnfColourPalette)
    return(sdsp)
}

##########################################
### Plotting boxplot (for mean length) ###
##########################################

#' Plotting performance: Box plot of mean width
#'
#' @param res The collected results from a computer experiment.
#' @return A plot object.
#' @export
lenBoxByType <- function(res) {
    lenp <- (  ggplot(res, aes(factor(confidence), width))
             + geom_boxplot(aes(fill=type))
             + coord_flip()
             + scale_x_discrete(name=expression("Aspired confidence level"))
             + scale_y_continuous(name=expression("Estimated mean interval width"))
             + scale_fill_discrete(expression("Type of method\nused for\ninterval estimation"))
             + theme(  legend.position="bottom"
                     , axis.text = element_text(colour = "black"))
             )
    return(lenp)
}

#' Plotting performance: Box plot of mean width
#'
#' @param res The collected results from a computer experiment.
#' @return A plot object.
#' @export
lenBoxByMethod <- function(res) {
    lenp <- (  ggplot(res, aes(factor(confidence), width))
             + geom_boxplot(aes(fill=method))
             + coord_flip()
             + scale_x_discrete(name=expression("Aspired confidence level"))
             + scale_y_continuous(name=expression("Estimated mean interval width"))
             + scale_fill_discrete(expression("Type of method\nused for\ninterval estimation"))
             + theme(  legend.position="bottom"
                     , axis.text = element_text(colour = "black"))
             )
    return(lenp)
}

############################################
### Plotting densities (for mean length) ###
############################################

#' Plotting performance: Density estimate of mean width
#'
#' By type.
#' @param res The collected results from a computer experiment.
#' @return A plot object.
#' @export
lenDenByType <- function(res) {
    lenp <- (  ggplot(res, aes(x=width))
             + geom_density(aes(fill=type), alpha=.3)
             + scale_x_continuous(name=expression("Mean interval width"))
             + scale_y_continuous(name=expression("Density"))
             + scale_fill_discrete(expression(paste("Type of ", tau, "-estimator")))
             + theme(  legend.position="bottom"
                     , axis.text = element_text(colour = "black"))
             )
    return(lenp)
}

#' Plotting performance: Density estimate of mean width
#'
#' By method.
#' @param res The collected results from a computer experiment.
#' @return A plot object.
#' @export
lenDenByMethod <- function(res) {
    lenp <- (  ggplot(res, aes(x=width))
             + geom_density(aes(fill=method), alpha=.3)
             + scale_x_continuous(name=expression("Mean interval width"))
             + scale_y_continuous(name=expression("Density"))
             + scale_fill_discrete(expression(paste("Method used for\ninterval estimation")))
             + theme(  legend.position="bottom"
                     , axis.text = element_text(colour = "black"))
             )
    return(lenp)
}
