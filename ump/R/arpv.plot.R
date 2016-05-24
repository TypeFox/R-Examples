arpv.plot <- function(alpha, phi, df = TRUE, verticals = TRUE) {

    if (! is.numeric(alpha)) stop("alpha not numeric")
    if (! is.numeric(phi)) stop("phi not numeric")
    if (! is.logical(df)) stop("df not logical")

    if (length(alpha) != length(phi)) stop("alpha and phi not same length")

    if (! all(0 <= alpha & alpha <= 1)) stop("alpha not in [0, 1]")
    if (! all(0 <= phi & phi <= 1)) stop("phi not in [0, 1]")

    if (df) {
        plot(alpha, phi, xlab = "significance level",
            ylab = "fuzzy P-value", type = "l")
        u <- par("usr")
        lines(c(u[1], alpha[1]), c(0, 0))
        lines(c(alpha[length(alpha)], u[2]), c(1, 1))
    } else {
        dens <- diff(phi) / diff(alpha)
        plot(range(alpha), range(0, dens), type = "n",
            xlab = "significance level",
            ylab = "density of randomized P-value")
        nalpha <- length(alpha)
        ndens <- length(dens)
        segments(alpha[-nalpha], dens, alpha[-1], dens)
        if (verticals) {
            # outer verticals
            jalpha <- c(1, nalpha)
            jdens <- c(1, ndens)
            segments(alpha[jalpha], rep(0, 2), alpha[jalpha], dens[jdens],
                lty = 2)
            # inner verticals
            if (nalpha > 2)
                segments(alpha[-jalpha], pmin(dens[-1], dens[-ndens]),
                    alpha[-jalpha], pmax(dens[-1], dens[-ndens]), lty = 2)
        }
    }
}
