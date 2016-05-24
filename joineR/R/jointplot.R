jointplot <- function (object, Y.col, Cens.col, lag, split = TRUE, col1, col2, xlab, ylab, gp1lab, gp2lab, smooth = 2/3, mean.profile = FALSE, mcol1, mcol2) 
{
    if (!is.vector(Y.col) | length(Y.col) > 1) {
        stop("Only one longitudinal response is possible to plot")
    }
    longdat <- object$longitudinal#[complete.cases(object$longitudinal), ]
    survdat <- object$survival
    if (is.numeric(Y.col)) {
        Y <- longdat[, Y.col]
    } else {
        Y <- longdat[[Y.col]]
        Y.col <- which(names(longdat) %in% Y.col)
    }
    t <- longdat[[object$time.col]]
    id <- longdat[[object$subj.col]]
    nobs <- diff(match(unique(id), id))
    nobs[length(nobs) + 1] <- length(id) - sum(nobs)
    index <- cumsum(nobs)
    cens <- survdat[, Cens.col]
    ft <- rep(t[index], nobs)
    t0 <- t - ft
    hue <- length(id)
    if (missing(lag)) {
        lag <- max(t)
    }
    if (missing(col1)) {
        col1 <- "blue"
    }
    if (missing(col2)) {
        col2 <- "red"
    }
    if (missing(xlab)) {
        xlab <- "Time"
    }
    if (missing(ylab)) {
        ylab <- "Y"
    }
    if (missing(gp1lab)) {
        gp1lab <- "Censored"
    }
    if (missing(gp2lab)) {
        gp2lab <- "Failed"
    }
    if (missing(mcol1)) {
        mcol1 <- "black"
    }
    if (missing(mcol2)) {
        mcol2 <- "black"
    }
    if (missing(smooth)) {
	  smooth <- 2/3
    }
    ii <- (cens == 0)
    hue[ii] <- col1
    hue[!ii] <- col2
    fac <- rep(cens, nobs)
    ii <- (fac == 0)
    fac[ii] <- gp1lab
    fac[!ii] <- gp2lab
    if (mean.profile == FALSE){
    if (split == TRUE) {
        xyplot(Y ~ t0 | fac, groups = id, type = "l", lty = 1, 
            xlim = c(-lag, 0), col = hue, xlab = xlab, ylab = ylab)
    }
    else {
        xyplot(Y ~ t0, groups = id, type = "l", xlim = c(-lag, 
            0), col = hue, xlab = xlab, ylab = ylab)
    }
} 
else {
    s1 <- lowess(t0[ii],Y[ii],f=smooth)
    s2 <- lowess(t0[!ii],Y[!ii],f=smooth)
    if(!is.character(summary(object)$times)){
	mean_cens <- unique(s1$y)
	mean_fail <-unique(s2$y)
	t_mean_cens <- unique(s1$x)
	t_mean_fail <- unique(s2$x)
    }
else {
    	mean_cens <- s1$y
	mean_fail <- s2$y
	t_mean_cens <- seq(min(t0[ii]),max(t0[ii]),length=length(mean_cens))
	t_mean_fail <- seq(min(t0[!ii]),max(t0[!ii]),length=length(mean_fail))
    }	

    Y <- c(Y,mean_cens,mean_fail)
    fac <- c(fac,rep(gp1lab,length(mean_cens)),rep(gp2lab,length(mean_fail)))
    id <- c(id,rep(max(id)+1,length(mean_cens)),rep(max(id)+2,length(mean_fail)))
    t0 <- c(t0,t_mean_cens,t_mean_fail)
    hue <- c(hue,mcol1,mcol2)

    if (split == TRUE) {
        xyplot(Y ~ t0 | fac, groups = id, type = "l", lty = 1, 
	lwd = c(rep(1,length(cens)),2,2),
        xlim = c(-lag, 0), col = hue, xlab = xlab, ylab = ylab, scales = list(alternating = FALSE))
    }
    else {
        xyplot(Y ~ t0, groups = id, type = "l", xlim = c(-lag, 
            0), col = hue, xlab = xlab, ylab = ylab)
    }
}
}
