predIntNormTestPower <-
function (n, df = n - 1, n.mean = 1, k = 1, delta.over.sigma = 0, 
    pi.type = "upper", conf.level = 0.95) 
{
    pi.type <- match.arg(pi.type, c("upper", "lower"))
    if (!is.vector(n, mode = "numeric") || !all(is.finite(n)) || 
        any(n < 2)) 
        stop(paste("'n' must be a numeric vector", "with all elements greater than or equal to 2", 
            "and no Missing (NA), Infinite(-Inf, Inf)", "or Undefined (Nan) values."))
    if (!is.vector(df, mode = "numeric") || !all(is.finite(df)) || 
        any(df < 1)) 
        stop(paste("'df' must be a numeric vector", "with all elements greater than or equal to 1", 
            "and no Missing (NA), Infinite(-Inf, Inf)", "or Undefined (Nan) values."))
    if (!is.vector(k, mode = "numeric") || !all(is.finite(k)) || 
        !all(k == trunc(k)) || any(k < 1) || !is.vector(n.mean, 
        mode = "numeric") || !all(is.finite(n.mean)) || !all(n.mean == 
        trunc(n.mean)) || any(n.mean < 1)) 
        stop(paste("'k' and 'n.mean' must be numeric vectors of integers", 
            "with all elements greater than 0", "and no Missing (NA), Infinite(-Inf, Inf),", 
            "or Undefined (Nan) values."))
    if (!is.vector(delta.over.sigma, mode = "numeric") || any(is.na(delta.over.sigma))) 
        stop(paste("'delta.over.sigma' must be a numeric vector", 
            "with no Missing (NA) or Undefined (Nan) values."))
    if (!is.vector(conf.level, mode = "numeric") || !all(is.finite(conf.level)) || 
        any(conf.level <= .Machine$double.eps) || any(conf.level >= 
        1 - .Machine$double.eps)) 
        stop(paste("'conf.level' must be a numeric vector", "with all elements between 0 and 1", 
            "and no Missing (NA), Infinite(-Inf, Inf),", "or Undefined (Nan) values."))
    arg.mat <- cbind.no.warn(n = as.vector(n), df = as.vector(df), 
        n.mean = as.vector(n.mean), k = as.vector(k), delta.over.sigma = as.vector(delta.over.sigma), 
        conf.level = as.vector(conf.level))
    for (i in c("n", "df", "n.mean", "k", "delta.over.sigma", 
        "conf.level")) assign(i, arg.mat[, i])
    N <- length(n)
    power <- numeric(N)
    index <- delta.over.sigma == 0
    if (any(index)) 
        power[index] <- 1 - conf.level[index]
    if (!all(index)) {
        for (i in c("n", "df", "n.mean", "k", "delta.over.sigma", 
            "conf.level")) assign(i, arg.mat[!index, i])
        if (all(k == 1) & all(n.mean == 1) & all(df == (n - 1))) {
            alternative <- ifelse(pi.type == "upper", "greater", 
                "less")
            power.sub <- tTestPower(n.or.n1 = n, n2 = 1, delta.over.sigma = delta.over.sigma, 
                alpha = 1 - conf.level, sample.type = "two.sample", 
                alternative = alternative, approx = FALSE)
        }
        else {
            N <- length(n)
            if (all(n == n[1]) & all(df == df[1]) & all(n.mean == 
                n.mean[1]) & all(k == k[1]) & all(conf.level == 
                conf.level[1])) {
                K <- predIntNormK(n = n[1], df = df[1], n.mean = n.mean[1], 
                  k = k[1], method = "exact", pi.type = pi.type, 
                  conf.level = conf.level[1])
                K <- rep(K, N)
            }
            else {
                K <- numeric(N)
                for (i in 1:N) {
                  K[i] <- predIntNormK(n = n[i], df = df[i], 
                    n.mean = n.mean[i], k = k[i], method = "exact", 
                    pi.type = pi.type, conf.level = conf.level[i])
                }
            }
            power.sub <- numeric(N)
            for (i in 1:N) {
                power.sub[i] <- predIntNormTestPowerScalar(n = n[i], 
                  df = df[i], K = K[i], n.mean = n.mean[i], k = k[i], 
                  delta.over.sigma = delta.over.sigma[i], pi.type = pi.type, 
                  conf.level = conf.level[i])
            }
        }
        power[!index] <- power.sub
    }
    names(power) <- NULL
    power
}
