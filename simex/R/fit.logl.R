fit.logl <-
function (lambda, p.names,estimates)
{
est <- (t(t(estimates) + (abs(apply(estimates, 2, min)) + 1) *
(apply(estimates, 2, min) <= 0)))
# calculating starting values
    start.mod <- lm(I(log(t(t(estimates) + (abs(apply(estimates, 2, min)) + 1) *
(apply(estimates, 2, min) <= 0)))) ~ lambda)
    start.val <- coef(start.mod)
    extrapolation <- list()
    # doing the extrapolation step
    for (d in p.names) {
        extrapolation[[d]] <- try(nls(est[, d] ~ exp(gamma.0 + (gamma.1 * lambda)),
start = list(gamma.0 = start.val[1, d], gamma.1 = start.val[2, d])),
silent = TRUE )
    }
    # security, in the case that nls() does not converge the "logl"-method is used
    if (any(lapply(extrapolation, class) == "try-error")) {
        warning("Function nls() did not converge, using 'logl' as fitting method",
call. =FALSE)
        extrapolation <- start.mod
    }
    return(extrapolation)
}

