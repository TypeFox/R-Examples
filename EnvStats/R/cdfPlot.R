cdfPlot <-
function (distribution = "norm", param.list = list(mean = 0, 
    sd = 1), left.tail.cutoff = ifelse(is.finite(supp.min), 0, 
    0.001), right.tail.cutoff = ifelse(is.finite(supp.max), 0, 
    0.001), plot.it = TRUE, add = FALSE, n.points = 1000, cdf.col = "black", 
    cdf.lwd = 3 * par("cex"), cdf.lty = 1, curve.fill = FALSE, 
    curve.fill.col = "cyan", digits = .Options$digits, ..., type = ifelse(discrete, 
        "s", "l"), main = NULL, xlab = NULL, ylab = NULL, xlim = NULL, 
    ylim = NULL) 
{
    if (!is.list(param.list)) 
        stop("'param.list' must be a list.")
    check.da.list <- check.distribution.args(distribution, param.list)
    dist.abb <- check.da.list$dist.abb
    dist.name <- check.da.list$dist.name
    dist.type <- check.da.list$dist.type
    n.dist.params <- check.da.list$n.dist.params
    dist.params.names <- check.da.list$dist.params.names
    param.list <- check.da.list$param.list
    supp.min <- eval(parse(text = EnvStats::Distribution.df[dist.abb, 
        "Support.Min"]), envir = param.list)
    supp.max <- eval(parse(text = EnvStats::Distribution.df[dist.abb, 
        "Support.Max"]), envir = param.list)
    qname <- paste("q", dist.abb, sep = "")
    pname <- paste("p", dist.abb, sep = "")
    if (left.tail.cutoff == 0) {
        if (supp.min == -Inf) 
            stop(paste("The value of 'left.tail.cutoff' must be greater", 
                "than 0 for the", dist.name, "distribution since the support", 
                "on the left-hand tail is infinite."))
        else q.min <- supp.min
    }
    else q.min <- do.call(qname, c(list(p = left.tail.cutoff), 
        param.list))
    if (right.tail.cutoff == 0) {
        if (supp.max == Inf) 
            stop(paste("The value of 'right.tail.cutoff' must be greater", 
                "than 0 for the", dist.name, "distribution since the support", 
                "on the right-hand tail is infinite."))
        else q.max <- supp.max
    }
    else q.max <- do.call(qname, c(list(p = 1 - right.tail.cutoff), 
        param.list))
    if (add) {
        usr <- par("usr")
        q.min <- max(usr[1], q.min)
        q.max <- min(usr[2], q.max)
    }
    if (any(dist.type == c("Continuous", "Mixed"))) {
        discrete <- FALSE
        q <- seq(q.min, q.max, len = n.points)
    }
    else {
        discrete <- TRUE
        q <- ceiling(q.min):floor(q.max)
    }
    y <- do.call(pname, c(list(q = q), param.list))
    if (plot.it) {
        if (!add) {
            gen.gp.list <- checkGraphicsPars(...)$gen.gp.list
            if (is.null(main)) {
                if (any(dist.abb == c("beta", "chisq", "f", "t"))) {
                  if (param.list$ncp > 0) 
                    main <- paste("Non-central ", dist.name, 
                      " CDF with\n", "(", paste(paste(dist.params.names, 
                        signif(unlist(param.list), digits), sep = "="), 
                        collapse = ", "), ")", sep = "")
                  else {
                    main <- paste(dist.name, " CDF with\n", "(", 
                      paste(paste(dist.params.names[-n.dist.params], 
                        signif(unlist(param.list[-n.dist.params]), 
                          digits), sep = "="), collapse = ", "), 
                      ")", sep = "")
                  }
                }
                else main <- paste(dist.name, " CDF with\n", 
                  "(", paste(paste(dist.params.names, signif(unlist(param.list), 
                    digits), sep = "="), collapse = ", "), ")", 
                  sep = "")
            }
            if (is.null(xlab)) 
                xlab <- "Value of Random Variable"
            if (is.null(ylab)) 
                ylab <- "Cumulative Frequency"
            if (is.null(xlim)) 
                xlim <- range(q)
            if (is.null(ylim)) 
                ylim <- c(0, 1)
            plot(q, y, type = "n", ..., xlab = xlab, ylab = ylab, 
                xlim = xlim, ylim = ylim, main = main)
            arg.list <- list(x = q, y = y)
            arg.list <- c(arg.list, gen.gp.list, list(type = type, 
                col = cdf.col, lwd = cdf.lwd, lty = cdf.lty))
            do.call("lines", arg.list)
        }
        else {
            lines(q, y, ..., type = type, col = cdf.col, lwd = cdf.lwd, 
                lty = cdf.lty)
        }
        if (discrete && left.tail.cutoff == 0) {
            lines(q[c(1, 1)], c(0, y[1]), ..., type = type, col = cdf.col, 
                lwd = cdf.lwd, lty = cdf.lty)
        }
        if (curve.fill) {
            n <- length(y)
            if (type != "s") {
                polygon(c(q, rev(q)), c(y, rep(0, n)), border = FALSE, 
                  col = curve.fill.col)
            }
            else {
                dum.q <- as.vector(matrix(q, 2, n, byrow = TRUE))
                m <- 2 * n
                dum.q <- dum.q[-c(1, m)]
                dum.y <- as.vector(matrix(y, 2, n, byrow = TRUE))
                dum.y <- dum.y[-c(m - 1, m)]
                polygon(c(dum.q, rev(dum.q)), c(dum.y, rep(0, 
                  m - 2)), border = FALSE, col = curve.fill.col)
            }
        }
    }
    invisible(list(Quantiles = q, Cumulative.Probabilities = y))
}
