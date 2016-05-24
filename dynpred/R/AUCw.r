#' Calculate dynamic AUC(t) curve
#' 
#' Calculate dynamic model-free curve of Area Under the Curve values over time,
#' based on the dynamic/incident AUC of Heagerty and Zheng.
#' 
#' 
#' @param formula Formula for prediction model to be used as in
#' \code{\link[survival:coxph]{coxph}}
#' @param data Data set in which to interpret the formula
#' @param width Width of the window
#' @return A data frame with columns \item{time}{The time points t at which
#' AUCw(t) changes value (either t or t+width is an event time point)}
#' \item{AUCw}{The AUCw(t) function} and with attribute \code{"width"} given as
#' input.
#' @author Hein Putter \email{H.Putter@@lumc.nl}
#' @references van Houwelingen HC, Putter H (2012). Dynamic Prediction in
#' Clinical Survival Analysis. Chapman & Hall.
#' @keywords univar
#' @examples
#' 
#' data(ova)
#' AUCw(Surv(tyears, d) ~ Karn + Broders + FIGO + Ascites + Diam, data = ova,
#'   width = 2)
#' 
#' @export AUCw
AUCw <- function(formula,data,width)
{
    cox1 <- coxph(formula, data)
    y <- cox1[["y"]]
    p <- ncol(y)
    time <- y[,p-1]
    status <- y[,p]
    x <- cox1$linear.predictors
    tt <- sort(unique(time[status==1]))
    ttw <- c(tt,tt-width)
    ttw <- ttw[ttw>0]
    ttw <- sort(unique(ttw))
    ntw <- length(ttw)
    AUCw <- rep(NA,ntw)
    for (j in 1:ntw) {
        twj <- ttw[j]
        ttj <- tt[((tt >= twj) & (tt <= twj+width))]
        ntj <- length(ttj)
        AUCt <- rep(NA,ntj)
        numsum <- denomsum <- 0
        for (i in 1:ntj) {
            ti <- ttj[i]
            # risk set
            Y <- sum(time>=ti)
            R <- which(time>ti) # !!! R(ti) is which(time>=ti), but for the "controls", j=i should be excluded, only valid without ties
            xi <- x[time==ti]
            num <- sum(x[R]<xi) + 0.5*sum(x[R]==xi)
            numsum <- numsum + num
            denomsum <- denomsum + Y-1 # Also only valid in absence of ties
        }
        AUCw[j] <- numsum/denomsum
    }
    res <- data.frame(time=ttw,AUCw=AUCw)
    attr(res, "width") <- width
    return(res)
}
