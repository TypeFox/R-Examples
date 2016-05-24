LCTS <-
function (cfs, tsx, tsy, filter.number = 1,
    family = c("DaubExPhase","DaubLeAsymm"), 
    plot.it = FALSE, spec.filter.number = 1, spec.family = c("DaubExPhase", "DaubLeAsymm")) 
{
family <- match.arg(family)
spec.family <- match.arg(spec.family)
#
# Compute time-varying linear combination of two time series and returns
# result of stationarity test on the combination, Z_t
#
#
# Compute dimensions of coefficients and time series and norm of coefficients
#
lcfs <- length(cfs)
lts <- length(tsx)
cfsnorm <- sum((cfs)^2)
#
# Separate out alpha and beta coefficients from single coefficient vector
# cfs and reform function version of both. 
#
alpha <- cfs[1:(lcfs/2)]
beta <- cfs[(lcfs/2 + 1):lcfs]
v <- coeftofn(alpha = alpha, beta = beta, n = lts,
	filter.number = filter.number, family = family)
#
# Form Z_t
#
lcts <- v$alpha * tsx + v$beta * tsy
#
# Compute EWS of lcts
#
lctsspec <- ewspec(lcts,
	filter.number = spec.filter.number, family = spec.family)$S
#
# Compute test statistic on spectrum
#
ans <- TOSts(lctsspec)/(cfsnorm^2)
#
# Plot the combination functions, Z_t and its spectrum if necessary.
#
if (plot.it == TRUE) {
	ts.plot(v$alpha, main = "alpha")
        ts.plot(v$beta, main = "beta")
        ts.plot(lcts, main = "Combined")
        plot(lctsspec, main = paste("Var is", signif(ans, 3)))
    }
#
# Return the test statistic
#
return(ans)
}
