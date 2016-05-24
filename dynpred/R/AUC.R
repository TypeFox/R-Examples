#' Calculate AUC(t) curve
#' 
#' Calculate model-free curve of Area Under the Curve values over time, based
#' on the dynamic/incident AUC of Heagerty and Zheng.
#' 
#' 
#' @param formula Formula for prediction model to be used as in
#' \code{\link[survival:coxph]{coxph}}
#' @param data Data set in which to interpret the formula
#' @param plot Determines whether the AUC function should be plotted (if
#' \code{TRUE} (default)) along with a \code{\link{lowess}} curve or not (if
#' \code{FALSE})
#' @return A list with elements \item{AUCt}{A data frame with time t in column
#' \code{time} and AUC(t) in column \code{AUC}} \item{AUC}{The AUC(t) weighted
#' by Y(t)-1, with Y(t) the number at risk at t; this coincides with Harrell's
#' c-index}
#' @author Hein Putter \email{H.Putter@@lumc.nl}
#' @references Harrell FE, Lee KL & Mark DB (1996), Multivariable prognostic
#' models: issues in developing models, evaluating assumptions and adequacy,
#' and measuring and reducing errors, Statistics in Medicine 15, 361-387.
#' 
#' Heagerty PJ & Zheng Y (2005), Survival model predictive accuracy and ROC
#' curves, Biometrics 61, 92-105.
#' 
#' van Houwelingen HC, Putter H (2012). Dynamic Prediction in Clinical Survival
#' Analysis. Chapman & Hall.
#' @keywords univar
#' @examples
#' 
#' data(ova)
#' AUC(Surv(tyears, d) ~ Karn + Broders + FIGO + Ascites + Diam, data = ova)
#' 
#' @export AUC
#' @import survival
#' @importFrom graphics lines abline
#' @importFrom stats lowess
AUC <- function(formula,data,plot=TRUE)
{
    cox1 <- coxph(formula, data)
    y <- cox1[["y"]]
    p <- ncol(y)
    time <- y[,p-1]
    status <- y[,p]
    x <- cox1$linear.predictors
    tt <- sort(unique(time[status==1]))
    nt <- length(tt)
    AUCt <- rep(NA,nt)
    numsum <- denomsum <- 0
    for (i in 1:nt) {
        ti <- tt[i]
        # risk set
        Y <- sum(time>=ti)
        R <- which(time>ti) # !!! R(ti) is which(time>=ti), but for the "controls", j=i should be excluded, only valid without ties
        xi <- x[time==ti]
        num <- sum(x[R]<xi) + 0.5*sum(x[R]==xi)
        AUCt[i] <- num/(Y-1) # Also only valid in absence of ties
        numsum <- numsum + num
        denomsum <- denomsum + Y-1 # Also only valid in absence of ties
    }
    AUC <- numsum/denomsum
    if (plot) {
        plot(tt,AUCt,xlab="Time t",ylab="AUC(t)")
        lines(lowess(data.frame(tt,AUCt)))
        abline(h=0.5,lty=3)
    }
    return(list(AUCt=data.frame(time=tt,AUC=AUCt),AUC=AUC))
}
