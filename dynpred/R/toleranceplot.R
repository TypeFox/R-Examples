#' Create a tolerance plot
#' 
#' Create a tolerance plot according to the methods of Henderson, Jones & Stare (2001)
#' 
#' Warnings will be issued each time the survival curve corresponding to a
#' value of x never goes below (1-coverage)/2; these warnings may be ignored.
#' 
#' @param formula Formula for prediction model to be used as in
#' \code{\link[survival:coxph]{coxph}}
#' @param data Data set in which to interpret the formula
#' @param coverage The coverage for the tolerance intervals (default is 0.8)
#' @param horizon The horizon, maximum value to be imputed in case of censored
#' observations; default is 1.05 times largest event time
#' @param plot Should the tolerance plot actually be plotted? Default is
#' \code{TRUE}
#' @param xlab Label for x-axis
#' @return A data frame with columns \item{x}{Predictor (centered at zero)}
#' \item{lower}{Lower bound of tolerance interval} \item{upper}{Upper bound of
#' tolerance interval} and with attributes \code{"coverage"} and
#' \code{"horizon"} (copied from input or default).
#' @author Hein Putter \email{H.Putter@@lumc.nl}
#' @references Henderson R, Jones M & Stare J (2001), Accuracy of point
#' predictions in survival analysis, Statistics in Medicine 20, 3083-3096.
#' 
#' van Houwelingen HC, Putter H (2012). Dynamic Prediction in Clinical Survival
#' Analysis. Chapman & Hall.
#' @keywords univar
#' @examples
#' 
#' data(ova)
#' toleranceplot(Surv(tyears, d) ~ Karn + Broders + FIGO + Ascites + Diam, data = ova)
#' 
#' @export toleranceplot
toleranceplot <- function(formula,data,coverage=0.8,horizon,plot=TRUE,xlab)
{
    cox1 <- coxph(formula, data)
    y <- cox1[["y"]]
    p <- ncol(y)
    time <- y[,p-1]
    status <- y[,p]
    x <- cox1$linear.predictors
    xs <- sort(unique(x))
    nx <- length(xs)
    if (missing(horizon)) horizon <- max(time[status==1])*1.05
    ylim <- c(0,horizon)
    if (missing(xlab)) xlab <- "x"
    if (plot) plot(range(xs),c(0,horizon),type="n",xlab=xlab,ylab="Tolerance interval")
    cx <- coxph(Surv(time,status) ~ x, method="breslow")
    res <- matrix(NA,nx,3)
    for (i in 1:nx) {
        xi <- xs[i]
        nd <- data.frame(x=xi)
        sf <- survfit(cx,newdata=nd)
        sf <- data.frame(time=sf$time,surv=sf$surv)
        low <- max(sf$time[sf$surv>1-(1-coverage)/2])
        up <- min(sf$time[sf$surv<(1-coverage)/2])
        if (is.infinite(up)) up <- horizon
        lines(rep(xi,2),c(low,up),type="l")
        res[i,] <- c(xi,low,up)
    }
    res <- as.data.frame(res)
    names(res) <- c("x","lower","upper")
    attr(res, "coverage") <- coverage
    attr(res, "horizon") <- horizon
    return(res)
}
