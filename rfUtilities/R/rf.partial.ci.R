#' @title Random Forests regression partial dependency plot with confidence intervals
#' @description Plots the partial dependency, and associated confidence intervals, of a random forests regression model 
#' 
#' @param m              randomForest regression object
#' @param x              data.frame or matrix of independent variables used to build model
#' @param yname          Name of the dependent variable
#' @param xname          Name of the independent variable for calculating partial dependence 
#' @param lci            Percentile of predictions to plot as the lower bound.
#' @param uci            Percentile of predictions to plot as the upper bound.
#' @param delta          Plot change in prediction between the independent variable and simulated values (Default = NULL)
#'
#' @return              recordedplot object to recall plot
#'
#' @note depends: randomForest
#'
#' @author Jeffrey S. Evans  <jeffrey_evans@@tnc.org>
#'                                                                   
#' @examples 
#'  library(randomForest)
#'  data(airquality)
#'  airquality <- na.omit(airquality)
#'  rf.ozone <- randomForest(y=airquality[,"Ozone"], airquality[,2:ncol(airquality)])
#'
#'  par(mfrow=c(2,2))
#'    for(i in c("Solar.R", "Wind", "Temp", "Day")){
#'      rf.partial.ci(m=rf.ozone, x=airquality, yname="Ozone", xname=i, delta=TRUE) 
#'    } 
#'
#' @export
rf.partial.ci <- function(m, x, yname, xname, lci=0.25, uci=0.75, delta=FALSE) {
  if (!inherits(m, "randomForest")) stop("m is not randomForest class object")
  if(m$type != "regression") stop("classification is not supported")
    conf.int <-(uci-lci)*100
    temp <- sort(x[, xname])
    y.hat.mean <- vector()
    y.hat.lb <- vector()
    y.hat.ub <- vector()
    y <- stats::predict(m, x)
  for (i in 1:length(temp)){
    x[, xname] <- temp[i]
    y.hat <- stats::predict(m, x)
    if (delta == TRUE){ y.hat <- y.hat - y }
    y.hat.mean[i] <- stats::weighted.mean(y.hat)
    y.hat.lb[i] <- stats::quantile(y.hat, lci)
    y.hat.ub[i] <- stats::quantile(y.hat, uci)
  }
  m.ci <- as.data.frame(cbind(temp, y.hat.mean, y.hat.lb, y.hat.ub))
  names(m.ci) <- c(xname, "y.hat.mean", "lci", "uci")
    y.lim=c(min(c(m.ci$uci, rev(m.ci$lci))), max(c(c(m.ci$uci, rev(m.ci$lci)))))  
    graphics::plot(m.ci[,xname], m.ci[,"y.hat.mean"], type = "n", ylim=y.lim, 
              xlab=xname, ylab=paste("Predicted values of", yname, sep=" ") )
      graphics::title(paste(paste("Partial Dependence of", yname, "on", xname),
         paste( paste(conf.int, "%", sep=""), "confidence interval"),sep="\n"))  
      graphics::polygon(c(m.ci[,xname], rev(m.ci[,xname])), c(m.ci$uci, rev(m.ci$lci)), col = "gray86")
      graphics::lines(m.ci[,xname], m.ci[,"y.hat.mean"], type = "b", pch = 20)
  return(Plot = grDevices::recordPlot())
} 
