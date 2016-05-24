ehyper <-
function (x, m = NULL, total = NULL, k, method = "mle") 
{
    if (!is.vector(x, mode = "numeric") || length(x) != 1 || 
        !is.finite(x) || x != trunc(x) || !is.vector(k, mode = "numeric") || 
        length(k) != 1 || !is.finite(k) || k != trunc(k)) 
        stop("Missing values not allowed for 'x' and 'k' and they must be integers")
    if ((is.null(m) & is.null(total)) || (!is.null(m) & !is.null(total))) 
        stop("You must supply 'm' or 'total' but not both")
    if (!is.null(m)) {
        if (!is.vector(m, mode = "numeric") || length(m) != 1 || 
            !is.finite(m) || m != trunc(m)) 
            stop("Missing values not allowed for 'm' and it must be an integer")
        est.m <- FALSE
    }
    else {
        if (!is.vector(total, mode = "numeric") || length(total) != 
            1 || !is.finite(total) || total != trunc(total)) 
            stop("Missing values not allowed for 'total' and it must be an integer")
        est.m <- TRUE
    }
    data.name <- deparse(substitute(x))
    method <- match.arg(method, c("mle", "mvue"))
    if (est.m) {
        if (x < 0 || x > k || total < 1 || k < 1 || k > total) 
            stop(paste("Illegal value(s) for 'x', 'total', and/or 'k'. ", 
                "'x' must be between 0 and 'k', 'total' must be at least 1,", 
                "and 'k' must be between 1 and 'total'."))
        m <- switch(method, mle = {
            m.hat <- ((total + 1) * x)/k
            ifelse(m.hat == trunc(m.hat), m.hat - 1, floor(m.hat))
        }, mvue = (total * x)/k)
    }
    else {
        if (x < 0 || x > min(k, m) || m < 0 || k < 1) 
            stop(paste("Illegal value(s) for 'x', 'm', and/or 'k'. ", 
                "'x' must be between 0 and min(k, m),", "'m' must be at least 0,", 
                "and 'k' must be at least 1."))
        if (method == "mvue") 
            stop("You must set method='mle' when estimating 'total'")
        total <- floor((k * m)/x)
    }
    dist.params <- c(m = m, n = total - m, k = k)
    method <- paste(method, "for", ifelse(est.m, "'m'", "'m+n'"))
    ret.list <- list(distribution = "Hypergeometric", sample.size = 1, 
        parameters = dist.params, n.param.est = 1, method = method, 
        data.name = data.name)
    oldClass(ret.list) <- "estimate"
    ret.list
}
