#' @export
#' @importFrom graphics par
#' @importFrom graphics plot
#' @importFrom graphics abline
#' @importFrom graphics lines
#' @importFrom graphics title
#' @importFrom utils modifyList
#' @importFrom stats qnorm

acfpacf<-function (x, nac, npac, datastr, ...)
{
    acfpacf_full <- function(x, nac, npac, plfg, acalpha, pacalpha,
        datastr, valcol, thrcol, thrmhcol) {
        nx = length(x)
        if (nac < 1 | npac < 1) {
            stop(" 'nac' or 'npac' must be positive")
        }
        if (nac > nx) {
            nac = nx
        }
        if (npac > nx) {
            npac = nx
        }
        ac <- acfpacf.acf(x, 3)
        ac = ac$acf
        ac = as.matrix(ac)
        pac <- acfpacf.pacf(x, npac)
        pac = pac$pacf
        pac = as.matrix(pac)
        thr = qnorm(1 - acalpha/2, 0, 1/sqrt(nx))
        thrmh = qnorm(1 - (acalpha/2)/nac, 0, 1/sqrt(nx))
        conf <- matrix(0, nac, 2)
        conf[, 1] = thr * matrix(1, nac)
        conf[, 2] = -conf[, 1]
        confmh <- matrix(0, nac, 2)
        confmh[, 1] = thrmh * matrix(1, nac)
        confmh[, 2] = -confmh[, 1]
        pthr = qnorm(1 - pacalpha/2, 0, 1/sqrt(nx))
        pthrmh = qnorm(1 - (pacalpha/2)/npac, 0, 1/sqrt(nx))
        pconf <- matrix(0, npac, 2)
        pconf[, 1] = pthr * matrix(1, npac)
        pconf[, 2] = -pconf[, 1]
        pconfmh <- matrix(0, npac, 2)
        pconfmh[, 1] = thrmh * matrix(1, npac)
        pconfmh[, 2] = -pconfmh[, 1]

   if (plfg) {
            par(mfrow = c(2, 1))
            ac = ac[seq(1, nac)]
            plot(seq(0, (nac - 1)), ac, xlab = "lags", ylab = "ACF",
                t = "h", lwd = 2, col = valcol, xlim = c(0, nac), ylim = c(0,1))
            abline(h = 0, col = "black")
            lines(seq(0, (nac - 1)), conf[, 1], col = thrcol)
            lines(seq(0, (nac - 1)), confmh[, 1], col = thrmhcol)
            lines(seq(0, (nac - 1)), conf[, 2], type = "l", col = thrcol)
            lines(seq(0, (nac - 1)), confmh[, 2], col = thrmhcol)
            title(paste("Usual ACF of", datastr, " for n=", nx, "alpha = ",
                acalpha))
            pac = pac[seq(1, npac)]
            plot(seq(0, (npac - 1)), pac, xlab = "n = no. samples between", ylab = "PACF",
                t = "h", lwd = 2, col = valcol, xlim = c(0, npac))
            abline(h = 0, col = "black")
            lines(seq(0, (npac - 1)), pconf[, 1], col = thrcol)
            lines(seq(0, (npac - 1)), pconfmh[, 1], col = thrmhcol)
            lines(seq(0, (npac - 1)), pconf[, 2], type = "l",
                col = thrcol)
            lines(seq(0, (npac - 1)), pconfmh[, 2], col = thrmhcol)
            title(paste("Usual PACF of", datastr, " for n=", nx, "alpha =  ",
                pacalpha))
        }
    }
    L <- modifyList(list(plfg = 1, acalpha = 0.05, pacalpha = 0.05,
        valcol = "red", thrcol = "green", thrmhcol = "blue"),
        list(x = x, nac = nac, npac = npac, datastr = datastr,
            ...))
    do.call(acfpacf_full, L)
}








