#' Calculate dynamic "death within window" curve
#' 
#' Calculate dynamic "death within window" curve, in other words, one minus
#' fixed width conditional survival curves, defined as P(T<=t+w|T>t), for a
#' fixed window width w.
#' 
#' "Die within window function" with window w, Fw(t) = P(T<=t+w|T>t), evaluated
#' at all time points t where the estimate changes value, and associated
#' pointwise confidence intervals (if \code{variance}=\code{TRUE}).
#' 
#' Both estimate and pointwise lower and upper confidence intervals are based
#' on the negative exponential of the Nelson-Aalen estimate of the cumulative
#' hazard, so P(T<=t+w|T>t) is estimated as exp(- int_t^t+w hatH_NA(s) ds),
#' with hatH_NA the non-parametric Nelson-Aalen estimate.
#' 
#' Note: in \code{object}, no event time points at or below zero allowed
#' 
#' @param object \code{\link[survival:survfit]{survfit}} object, use
#' type="aalen"
#' @param width Width of the window
#' @param variance Boolean (default=\code{TRUE}); should pointwise confidence
#' interval of the probabilities be calculated?
#' @param conf.level The confidence level, between 0 and 1 (default=0.95)
#' @return A data frame with columns \item{time}{The time points t at which
#' Fw(t) changes value (either t or t+width is an event time point)}
#' \item{Fw}{The Fw(t) function} \item{low}{Lower end of confidence interval}
#' \item{up}{Upper end of confidence interval} and with attribute
#' \code{"width"} as given as input.
#' @author Hein Putter \email{H.Putter@@lumc.nl}
#' @references van Houwelingen HC, Putter H (2012). Dynamic Prediction in
#' Clinical Survival Analysis. Chapman & Hall.
#' @keywords univar
#' @examples
#' 
#' data(wbc1)
#' c0 <- coxph(Surv(tyears, d) ~ 1, data = wbc1, method="breslow")
#' sf0 <- survfit(c0)
#' Fw <- Fwindow(sf0,4)
#' 
#' @importFrom stats qnorm
#' 
#' @export Fwindow
Fwindow <- function(object, width, variance=TRUE, conf.level=0.95)
{
    if (variance)
        sf <- data.frame(time=object$time,surv=object$surv,varHaz=object$std.err^2)
    else sf <- data.frame(time=object$time,surv=object$surv)
    sf$Haz <- -log(sf$surv)
    tt <- c(0,sf$time) # assuming no event at t=0 or t<0
    ttt <- c(0,tt,tt-width)
    ttt <- ttt[ttt >= 0]
    ttt <- sort(unique(ttt))
    ttt <- unique(ttt)
    H <- outer(c(0,sf$Haz),c(0,sf$Haz),"-")
    dimnames(H) <- list(tt,tt)
    tt <- c(tt,Inf)
    idx1 <- as.numeric(cut(ttt,tt,right=FALSE))
    idx2 <- as.numeric(cut(ttt+width,tt,right=FALSE))
    Fw <- diag(H[idx2,idx1])
    nt <- length(Fw)
    Fw[nt] <- Fw[nt-1]
    if (variance) {
        varH <- outer(c(0,sf$varHaz),c(0,sf$varHaz),"-")
        varFw <- diag(varH[idx2,idx1])
        varFw[nt] <- varFw[nt-1]
        ciwidth <- qnorm(1-(1-conf.level)/2)*sqrt(varFw)
        low <- Fw - ciwidth
        up <- Fw + ciwidth
        low[low<0] <- 0 # negative lower values of cum hazard set to zero
    }
    # Return on probability scale
    Fw <- 1 - exp(-Fw)
    if (variance) {
        low <- 1 - exp(-low)
        up <- 1 - exp(-up)
        res <- data.frame(time=ttt,Fw=Fw,low=low,up=up)
    }
    else res <- data.frame(time=ttt,Fw=Fw)
    attr(res,"width") <- width
    return(res)
}
