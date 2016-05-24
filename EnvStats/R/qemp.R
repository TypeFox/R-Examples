qemp <-
function (p, obs, discrete = FALSE, prob.method = ifelse(discrete, 
    "emp.probs", "plot.pos"), plot.pos.con = 0.375) 
{
    if (!is.vector(p, mode = "numeric") || (length.p <- length(p)) < 
        1) 
        stop("'p' must be a numeric vector with at least one element")
    names.p <- names(p)
    if (!is.vector(obs, mode = "numeric")) 
        stop("'obs' must be a numeric vector")
    if ((bad.obs <- sum(!(obs.ok <- is.finite(obs)))) > 0) {
        is.not.finite.warning(obs)
        obs <- obs[obs.ok]
        warning(paste(bad.obs, "observations with NA/NaN/Inf in 'obs' removed."))
    }
    if (length(obs) < 2) 
        stop("'obs' must be a numeric vector with at least 2 non-missing elements")
    if (!is.vector(discrete, mode = "logical") || length(discrete) != 
        1) 
        stop("'discrete' must be a logical scalar")
    prob.method <- match.arg(prob.method, c("emp.probs", "plot.pos"))
    if (!is.vector(plot.pos.con, mode = "numeric") || length(plot.pos.con) != 
        1 || plot.pos.con < 0 || plot.pos.con > 1) 
        stop("'plot.pos.con' must be a numeric scalar between 0 and 1")
    na.index <- is.na(p)
    if (all(na.index)) 
        q <- as.numeric(rep(NA, length.p))
    else {
        q <- numeric(length.p)
        q[na.index] <- NA
        q.no.na <- q[!na.index]
        p.no.na <- p[!na.index]
        if (any(p.no.na < 0) || any(p.no.na > 1)) 
            stop("All non-missing values of 'p' must be between 0 and 1.")
        ecdf.list <- ecdfPlot(obs, discrete = discrete, prob.method = prob.method, 
            plot.pos.con = plot.pos.con, plot.it = FALSE)
        x <- ecdf.list$Order.Statistics
        F.x <- ecdf.list$Cumulative.Probabilities
        if (discrete) {
            q.no.na <- approx(x = F.x, y = x, xout = p.no.na, 
                method = "constant", f = 1, rule = 2)$y
        }
        else {
            q.no.na <- approx(x = F.x, y = x, xout = p.no.na, 
                method = "linear", rule = 2)$y
        }
        q[!na.index] <- q.no.na
    }
    if (!is.null(names.p)) 
        names(q) <- rep(names.p, length = length(q))
    else names(q) <- NULL
    q
}
