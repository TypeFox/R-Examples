#' Create scatter plot with imputed survival times
#' 
#' Create scatter plot with imputed survival times.
#' 
#' Imputation is used for censored survival times.
#' 
#' @param formula Formula for prediction model to be used as in
#' \code{\link[survival:coxph]{coxph}}
#' @param data Data set in which to interpret the formula
#' @param horizon The horizon, maximum value to be imputed in case of censored
#' observations; default is 1.05 times largest event time
#' @param plot Should the tolerance plot actually be plotted? Default is
#' \code{TRUE}
#' @param xlab Label for x-axis
#' @return A data frame with columns \item{x}{Predictor (centered at zero)}
#' \item{imputed}{(Imputed) survival time} and with attribute \code{"horizon"}
#' (copied from input or default).
#' @author Hein Putter \email{H.Putter@@lumc.nl}
#' @references Royston P (2001), The lognormal distribution as a model for
#' survival time in cancer, with an emphasis on prognostic factors, Statistica
#' Neerlandica 55, 89-104.
#' 
#' van Houwelingen HC, Putter H (2012). Dynamic Prediction in Clinical Survival
#' Analysis. Chapman & Hall.
#' @keywords univar
#' @examples
#' 
#' data(ova)
#' scatterplot(Surv(tyears, d) ~ Karn + Broders + FIGO + Ascites + Diam, data = ova)
#' 
#' @importFrom stats runif lsfit lm 
#' @importFrom graphics points text
#' 
#' 
#' @export scatterplot
scatterplot <- function(formula,data,horizon,plot=TRUE,xlab)
{
    cox1 <- coxph(formula, data)
    y <- cox1[["y"]]
    p <- ncol(y)
    time <- y[,p-1]
    status <- y[,p]
    x <- cox1$linear.predictors
    if (missing(horizon)) horizon <- max(time[status==1])*1.05
    ylim <- c(0,1.05*horizon) # extra space to print R-square
    imputed <- time
    cens <- which(status==0)
    events <- which(status==1)
    ncens <- length(cens)
    if (missing(xlab)) xlab <- "x"
    cx <- coxph(Surv(time,status) ~ x, method="breslow")
    for (i in 1:ncens) {
        xi <- x[cens[i]]
        nd <- data.frame(x=xi)
        sf <- survfit(cx,newdata=nd)
        sf <- data.frame(time=sf$time,surv=sf$surv)
        sf <- sf[!duplicated(sf$surv),]
        if (imputed[cens[i]]>max(sf$time))
            imputed[cens[i]] <- horizon
        if (imputed[cens[i]]<max(sf$time)) {
            sf1 <- sf[sf$time <= imputed[cens[i]],]
            sf2 <- sf[sf$time >= imputed[cens[i]],]
            sf2$surv <- 1-sf2$surv/min(sf1$surv)
            if (sf2$surv[1]==0) sf2 <- sf2[-1,]
            rand <- runif(1)
            survtimes <- c(0,sf2$surv,1)
            idx <- as.numeric(cut(rand,survtimes))
            if (idx>nrow(sf2)) imputed[cens[i]] <- horizon
            if (idx<=nrow(sf2)) imputed[cens[i]] <- sf2$time[idx]
        }
    }
    res <- data.frame(x=x, imputed=imputed)
    attr(res, "horizon") <- horizon
    if (plot) {
        plot(x, imputed, type="n", ylim=ylim, xlab=xlab, ylab="(Imputed) survival time")
        points(x[cens], imputed[cens])
        points(x[events], imputed[events], pch=16)
        abline(lsfit(x,imputed),lty=2)
        lmova <- lm(imputed ~ x)
        r2 <- summary(lmova)$r.squared
        text(max(x),1.05*horizon,paste("R-squared =",round(r2,3)),adj=1)
    }
    return(res)
}
