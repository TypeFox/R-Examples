#' Calculate prediction error curve
#' 
#' Calculate prediction error curve.
#' 
#' The \code{censformula} is used to calculate inverse probability of censoring
#' weights (IPCW).
#' 
#' @aliases pe pecox
#' @param time Vector of time points in data
#' @param status Vector of event indicators in data
#' @param tsurv Vector of time points corresponding to the estimated survival
#' probabilities in \code{survmat}
#' @param survmat Matrix of estimated survival probabilities; dimension should
#' be length of tsurv x length of time
#' @param tcens Vector of time points corresponding to the estimated censoring
#' probabilities in \code{censmat}
#' @param censmat Matrix of estimated censoring probabilities; dimension should
#' be length of tcens x length of time
#' @param FUN The error function, either \code{"KL"} (default) for
#' Kullback-Leibler or \code{"Brier"} for Brier score
#' @param tout Vector of time points at which to evaluate prediction error. If
#' missing, prediction error will be evaluated at all time points where the
#' estimate will change value
#' @param formula Formula for prediction model to be used as in
#' \code{\link[survival:coxph]{coxph}}
#' @param censformula Formula for censoring model, also to be used as in
#' \code{\link[survival:coxph]{coxph}}
#' @param data Data set in which to interpret \code{formula}
#' @param censdata Data set in which to interpret \code{censformula}
#' @param CV Boolean (default=\code{FALSE}); if \code{TRUE}, (leave-one-out)
#' cross-validation is used for the survival probabilities
#' @param progress Boolean (default=\code{FALSE}); if \code{TRUE}, progress is
#' printed on screen
#' @return A data frame with columns \item{time}{Event time points}
#' \item{Err}{Prediction error of model specified by \code{formula} at these
#' time points}
#' @author Hein Putter \email{H.Putter@@lumc.nl}
#' @references Graf E, Schmoor C, Sauerbrei W & Schumacher M (1999), Assessment
#' and comparison of prognostic classification schemes for survival data,
#' Statistics in Medicine 18, 2529-2545.
#' 
#' Gerds & Schumacher (2006), Consistent estimation of the expected Brier score
#' in general survival models with right-censored event times, Biometrical
#' Journal 48, 1029-1040.
#' 
#' van Houwelingen HC, Putter H (2012). Dynamic Prediction in Clinical Survival
#' Analysis. Chapman & Hall.
#' @keywords univar
#' @examples
#' 
#' data(ova)
#' # Example on a subset, because the effect of CV is clearer
#' ova2 <- ova[1:100,]
#' pecox(Surv(tyears, d) ~ Karn + Broders + FIGO + Ascites + Diam, Surv(tyears, 1-d) ~ 1,
#'   data = ova2, FUN="Brier", tout=seq(0,6,by=0.5))
#' pecox(Surv(tyears, d) ~ Karn + Broders + FIGO + Ascites + Diam, Surv(tyears, 1-d) ~ 1,
#'   data = ova2, FUN="Brier", tout=seq(0,6,by=0.5), CV=TRUE, progress=TRUE)
#' 
#' \donttest{
#' pecox(Surv(tyears, d) ~ Karn + Broders + FIGO + Ascites + Diam, Surv(tyears, 1-d) ~ 1,
#'   data = ova, FUN="Brier", tout=seq(0,6,by=0.5))
#' pecox(Surv(tyears, d) ~ Karn + Broders + FIGO + Ascites + Diam, Surv(tyears, 1-d) ~ 1,
#'   data = ova, FUN="Brier", tout=seq(0,6,by=0.5), CV=TRUE, progress=TRUE)
#' }
#' 
#' @export pe
pe <- function(time,status,tsurv,survmat,tcens,censmat,FUN=c("KL","Brier"),tout)
{
    FUN <- match.arg(FUN)
    FUNn <- 2
    if (FUN=="Brier") FUNn <- 1
    n <- length(time) # no checks that this coincides with ncols of survmat and censmat
    nsurv <- length(tsurv) # no check that this coincides with nrows of survmat
    ncens <- length(tcens) # no check that this coincides with nrows of censmat
    nout <- length(tout)
    res <- .C('prederr',
      sn = as.integer(n),
      time = as.double(time),
      status = as.integer(status),
      snsurv = as.integer(nsurv),
      sncens = as.integer(ncens),
      snout = as.integer(nout),
      tsurv = as.double(tsurv),
      survmat = as.double(survmat),
      tcens = as.double(tcens),
      censmat = as.double(censmat),
      tout = as.double(tout),
      sFUN = as.integer(FUNn),
      err = double(nout),
      work = double(n)
    )
    res <- data.frame(time=tout,Err=res$err)
    attr(res,"score") <- FUN
    return(res)
}

#' @rdname pe
#' @export
#' @useDynLib dynpred

pecox <- function(formula, censformula, data, censdata,
    FUN = c("KL", "Brier"), tout, CV = FALSE, progress = FALSE)
{
    cox1 <- coxph(formula, data)
    y <- cox1[["y"]]
    p <- ncol(y)
    time <- y[, p - 1]
    status <- y[, p]
    n <- length(time)
    if (nrow(data) != n)
        stop("missing data in time or status not allowed")
    ord <- order(time, -status)
    time <- time[ord]
    status <- status[ord]
    tt <- sort(unique(time[status == 1]))
    nt <- length(tt)
    if (!CV) {
        progress <- FALSE
        x <- cox1$linear.predictors[ord]
        cox1 <- coxph(Surv(time, status) ~ x)
        # if no covariates, then use Kalbfleisch-Prentice type in survfit
        if (sum(x^2)==0)
          sf <- survfit(cox1, newdata = data.frame(x = x), type="kalbfl")
        else
          sf <- survfit(cox1, newdata = data.frame(x = x))
        tt <- sf$time
        survmat <- sf$surv
    }
    if (CV) {
        x <- cox1$linear.predictors[ord]
        if (sum(x^2)==0) stop("No cross-validation for null model implemented")
        X <- model.matrix(formula, data = data)
        X <- X[, -1, drop = FALSE]
        X <- X[ord, , drop = FALSE]
        if (progress) {
            m <- floor(log10(n)) + 1
            pre <- rep("\b", 2 * m + 1)
            cat("Calculating cross-validated survival curves:\n")
        }
        survmat <- matrix(NA, nt, n)
        for (i in 1:n) {
            if (progress) {
                cat(pre, i, "/", n, sep = "")
                flush.console()
            }
            cmini <- coxph(Surv(time[-i], status[-i]) ~ X[-i,
                , drop = FALSE], method = "breslow")
            xi <- as.vector(X[-i, , drop = FALSE] %*% cmini$coef)
            ximini <- as.numeric(X[i, ] %*% cmini$coef)
            cmini <- coxph(Surv(time[-i], status[-i]) ~ xi, method = "breslow")
            survi <- survfit(cmini, newdata = data.frame(xi = ximini))
            survmat[, i] <- evalstep(survi$time, survi$surv,
                tt, subst = 1)
        }
        if (progress) {
            cat("\nCalculating prediction error ...")
            flush.console()
        }
    }
    if (tt[1] > 0) {
        tsurv <- c(0, tt)
        survmat <- rbind(rep(1, n), survmat)
    }
    else tsurv <- tt
    nsurv <- length(tsurv)
    ## censoring
    if (missing(censdata)) censdata <- data
    coxcens <- coxph(censformula, censdata)
    ycens <- coxcens[["y"]]
    p <- ncol(ycens)
    tcens <- ycens[, p - 1]
    dcens <- ycens[, p]
    xcens <- coxcens$linear.predictors
    coxcens <- coxph(Surv(tcens, dcens) ~ xcens)
    # if no covariates, then use Kalbfleisch-Prentice type in survfit
    if (sum(xcens^2)==0)
      sfcens <- survfit(coxcens, newdata = data.frame(xcens = xcens),
        type="kalbfl")
    else
      sfcens <- survfit(coxcens, newdata = data.frame(xcens = xcens))
    tcens <- sfcens$time
    censmat <- sfcens$surv
    if (tcens[1] > 0) {
        tcens <- c(0, tcens)
        censmat <- rbind(rep(1, n), censmat)
    }
    ncens <- length(tcens)
    if (missing(tout))
        tout <- unique(c(tsurv, tcens))
    tout <- sort(tout)
    nout <- length(tout)
    res <- pe(time, status, tsurv, survmat, tcens, censmat, FUN,
        tout)
    if (progress)
        cat("\n")
    return(res)
}
