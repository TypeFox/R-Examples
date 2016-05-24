ADCFplot <- function (x, MaxLag = 15, ylim=NULL, main = NULL, method = c("Wild Bootstrap", 
    "Subsampling"), b = 499) 
{
    if (b <= 0) 
        stop("No plot is given for b<=0")
    if (MaxLag==0)
        stop("MaxLag must be greater than 0")
    series <- deparse(substitute(x))
    method <- match.arg(method)
    if (missing(method)) 
        method = "Wild Bootstrap"
    n <- length(x)
    adcor <- ADCF(x, MaxLag)
    if (method == "Wild Bootstrap") {
        cv <- RbootCV(n, MaxLag, b = b, parallel = TRUE)
    }
    else {
        if (((n - MaxLag) < 0) || ((n - MaxLag) < 4) || ((n - 
            MaxLag) <= 25)) 
            stop("Give bigger sample size n")
        cv <- SubsCV(x, MaxLag, parallel = TRUE)
    }
    r1 <- max(cv, 1)
    if (is.null(ylim)) ylim=c(0,r1)
    if (length(cv) == 1) {
        if (is.null(main)) {
            plot(0:MaxLag, adcor, type = "n", main = paste("Series", 
                series), xlab = "Lag", ylab = "ADCF", ylim = ylim)
        }
        else {
            plot(0:MaxLag, adcor, type = "n", main = main, xlab = "Lag", 
                ylab = "ADCF", ylim = ylim)
        }
        for (i in seq(0, MaxLag, by = 1)) {
            segments(i, 0, i, adcor[i + 1])
        }
        points(0:MaxLag, rep(cv, MaxLag + 1), type = "l", lty = 3, 
            lwd = 2, col = "blue")
    }
    else {
        if (is.null(main)) {
            plot(1:MaxLag, adcor[-1], type = "n", main = paste("Series", 
                series), xlab = "Lag", ylab = "ADCF", ylim = ylim)
        }
        else {
            plot(1:MaxLag, adcor[-1], type = "n", main = main, 
                xlab = "Lag", ylab = "ADCF", ylim = ylim)
        }
        for (i in seq(1, MaxLag, by = 1)) {
            segments(i, 0, i, adcor[i + 1])
        }
        points(1:MaxLag, cv, type = "l", lty = 3, lwd = 2, col = "blue")
    }
    result <- list(ADCF = adcor, method = method, critical.values = cv)
    return(result)
}