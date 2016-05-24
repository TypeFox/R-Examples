pemp <-
function (q, obs, discrete = FALSE, prob.method = ifelse(discrete, 
    "emp.probs", "plot.pos"), plot.pos.con = 0.375) 
{
    if (!is.vector(q, mode = "numeric") || (length.q <- length(q)) < 
        1) 
        stop("'q' must be a numeric vector with at least one element")
    names.q <- names(q)
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
    na.index <- is.na(q)
    if (all(na.index)) 
        p <- as.numeric(rep(NA, length.q))
    else {
        p <- numeric(length.q)
        p[na.index] <- NA
        p.no.na <- p[!na.index]
        q.no.na <- q[!na.index]
        ecdf.list <- ecdfPlot(obs, discrete = discrete, prob.method = prob.method, 
            plot.pos.con = plot.pos.con, plot.it = FALSE)
        x <- ecdf.list$Order.Statistics
        F.x <- ecdf.list$Cumulative.Probabilities
        if (discrete) {
            min.x <- min(x)
            max.x <- max(x)
            p.no.na[q.no.na < min.x] <- 0
            p.no.na[q.no.na >= max.x] <- 1
            if (any(index <- min.x <= q.no.na & q.no.na < max.x)) {
                qs <- q.no.na[index]
                x.index <- 1:length(x)
                q.index <- sapply(1:length(qs), function(i, x.index, 
                  qs, x) max(x.index[qs[i] >= x]), x.index = x.index, 
                  qs = qs, x = x)
                p.no.na[index] <- F.x[q.index]
            }
        }
        else {
            p.no.na <- approx(x = x, y = F.x, xout = q.no.na, 
                method = "linear", rule = 2)$y
        }
        p[!na.index] <- p.no.na
    }
    if (!is.null(names.q)) 
        names(p) <- rep(names.q, length = length(p))
    else names(p) <- NULL
    p
}
