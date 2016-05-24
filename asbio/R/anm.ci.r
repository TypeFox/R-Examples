anm.ci <- function (parent=expression(rnorm(n)), par.val, conf = 0.95, sigma = NULL, par.type = c("mu", 
    "median", "sigma.sq", "p"), n.est = 100, n = 50, err.col = 2, par.col = 4, interval = 0.1, 
    ...) 
{
    ci <- matrix(ncol = 3, nrow = n.est)
    names <- c(expression(mu), "Pop. Median", expression(sigma^2), 
       expression(pi))
    for (i in 1:n.est) {
        x <- sample(eval(parent), size = n, replace = FALSE)
        if (par.type == "mu") {
            if (!is.null(sigma)) {
                cint <- ci.mu.z(x, conf, sigma)
                ci[i, ] <- c(cint$ci[1], cint$ci[2], cint$ci[3])
            }
            if (is.null(sigma)) {
                cint <- ci.mu.t(x, conf)
                ci[i, ] <- c(cint$ci[1], cint$ci[2], cint$ci[3])
            }
        }
        if (par.type == "median") {
            cint <- ci.median(x, conf)
            ci[i, ] <- c(cint$ci[1], cint$ci[2], cint$ci[3])
        }
        if (par.type == "sigma.sq") {
            cint <- ci.sigma(x, conf)
            ci[i, ] <- c(cint$ci[1], cint$ci[2], cint$ci[3])
        }
        if (par.type == "p") {
            cint <- ci.p(x, conf)
            ci[i, ] <- c(cint$ci[1], cint$ci[2], cint$ci[3])
        }
    }
    lcol <- matrix(nrow = n.est, ncol = 1)
    for (i in 1:n.est) {
        lcol[i] <- ifelse(ci[, 2][i] < par.val & ci[, 3][i] > 
            par.val, 1, err.col)
    }
    dev.hold()
    plot(ci[, 1], seq(1, n.est), xlim = c(min(ci[, 2]), max(ci[, 
        3])), type = "n", xlab = "Point and interval estimates", 
        ylab = "Estimate number", ...)
    abline(v = par.val, lty = 2, col = par.col)
    if (par.type == "mu") 
        mtext(names[1], 3, at = par.val, font = 3, line = 0.1)
    if (par.type == "median") 
        mtext(names[2], 3, at = par.val, font = 3, line = 0.1)
    if (par.type == "sigma.sq") 
        mtext(names[3], 3, at = par.val, font = 3, line = 0.1)
    if (par.type == "p") 
        mtext(names[4], 3, at = par.val, font = 3, line = 0.1)
    for (i in 1:n.est) {
        points(ci[, 1][i], i, pch = 19, cex = 0.6, col = lcol[i])
        segments(x0 = ci[, 2][i], x1 = ci[, 3][i], y0 = i, y1 = i, 
            col = lcol[i])
        dev.flush()
        Sys.sleep(interval)
    }
    mtext(bquote(paste("Obs. cvg. = ", .(round(sum(sapply(lcol == 
        1, sum))/n.est, 2)))), 3, at = max(ci[, 3]), adj = 1, 
        line = 1)
    mtext(bquote(paste("Conf. = ", .(conf))), 3, at = min(ci[, 
        2]), adj = 0, line = 1)
}
