predIntNormTestPowerScalar <-
function (n, df = n - 1, n.mean = 1, K, k = 1, delta.over.sigma = 0, 
    pi.type = c("upper", "lower"), conf.level = 0.95) 
{
    if (!is.vector(K, mode = "numeric") || length(K) != 1 || 
        K <= 0) 
        stop("'K' must be a positive scalar")
    if (!is.vector(n, mode = "numeric") || length(n) != 1 || 
        !is.vector(df, mode = "numeric") || length(df) != 1 || 
        !is.vector(k, mode = "numeric") || length(k) != 1 || 
        !is.vector(delta.over.sigma, mode = "numeric") || length(delta.over.sigma) != 
        1 || !is.vector(conf.level, mode = "numeric") || length(conf.level) != 
        1) 
        stop(paste("'n', 'df', 'k', 'delta.over.sigma',", "and 'conf.level' must be numeric scalars"))
    if (!is.finite(n) || n < 2) 
        stop("'n' must be finite and greater than 1")
    if (!is.finite(df) || df < 1) 
        stop("'df' must be finite and greater than 0")
    if (!is.finite(k) || k < 1 || k != trunc(k)) 
        stop("'k' must be an integer greater than or equal to 1")
    if (!is.finite(n.mean) || n.mean < 1 || n.mean != trunc(n.mean)) 
        stop("'n.mean' must be an integer greater than or equal to 1")
    if (is.na(delta.over.sigma)) 
        stop(paste("'delta.over.sigma' cannot be a missing (NA)", 
            "or undefined (NaN) value"))
    if (conf.level <= 0 || conf.level >= 1) 
        stop("'conf.level' must be between 0 and 1")
    pi.type <- match.arg(pi.type)
    if (delta.over.sigma == 0) 
        power <- 1 - conf.level
    else if (!is.finite(delta.over.sigma)) {
        if (pi.type == "upper") {
            if (delta.over.sigma == Inf) 
                power <- 1
            else power <- 0
        }
        else {
            if (delta.over.sigma == -Inf) 
                power <- 1
            else power <- 0
        }
    }
    else {
        if (k == 1) {
            dum <- sqrt((n.mean * n)/(n.mean + n))
            ncp <- dum * delta.over.sigma
            power <- switch(pi.type, upper = pT(q = dum * K, 
                df = df, ncp = ncp, lower.tail = FALSE), lower = pT(q = -dum * 
                K, df = df, ncp = ncp, lower.tail = FALSE))
        }
        else {
            power <- switch(pi.type, upper = 1 - predIntNormSimultaneousProb(n = n, 
                df = df, n.mean = n.mean, K = K, delta.over.sigma = delta.over.sigma, 
                k = k, m = k, r = 1, rule = "k.of.m"), lower = 1 - 
                predIntNormSimultaneousProb(n = n, df = df, n.mean = n.mean, 
                  K = K, delta.over.sigma = -delta.over.sigma, 
                  k = k, m = k, r = 1, rule = "k.of.m"))
        }
    }
    names(power) <- NULL
    power
}
