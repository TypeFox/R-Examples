sraFormatlegend <-
function (names, values, AIC = NULL, ...) 
{
    ans <- NULL
    for (nn in names) {
        if (nn == "mu0") {
            ans <- c(ans, substitute(mu[0] ~ ~vv, list(vv = format(values[nn], 
                ...))))
        }
        if (nn == "logvarA0") {
            ans <- c(ans, substitute({
                sigma^2
            }[A[0]] ~ ~vv, list(vv = format(exp(values[nn]), 
                ...))))
        }
        if (nn == "logIA0") {
            ans <- c(ans, substitute(I[A[0]] ~ ~vv, list(vv = format(exp(values[nn]), 
                ...))))
        }
        if (nn == "logith20") {
            ans <- c(ans, substitute({
                h^2
            }[0] ~ ~vv, list(vv = format(exp(values[nn])/(1 + 
                exp(values[nn])), ...))))
        }
        if (nn == "logvarE0") {
            ans <- c(ans, substitute({
                sigma^2
            }[E[0]] ~ ~vv, list(vv = format(exp(values[nn]), 
                ...))))
        }
        if (nn == "logvarME") {
            ans <- c(ans, substitute({
                sigma^2
            }[me] ~ ~vv, list(vv = format(exp(values[nn]), ...))))
        }
        if (nn == "logIE0") {
            ans <- c(ans, substitute(I[E[0]] ~ ~vv, list(vv = format(exp(values[nn]), 
                ...))))
        }
        if (nn == "logvarP0") {
            ans <- c(ans, substitute({
                sigma^2
            }[P[0]] ~ ~vv, list(vv = format(exp(values[nn]), 
                ...))))
        }
        if (nn == "logepsilon") {
            ans <- c(ans, substitute(epsilon ~ ~vv, list(vv = format(exp(values[nn]), 
                ...))))
        }
        if (nn == "logminusepsilon") {
            ans <- c(ans, substitute(epsilon ~ ~vv, list(vv = format(-exp(values[nn]), 
                ...))))
        }
        if (nn == "logvarepsilon") {
            ans <- c(ans, substitute({
                sigma^2
            }[epsilon] ~ ~vv, list(vv = format(exp(values[nn]), 
                ...))))
        }
        if (nn == "kc") {
            ans <- c(ans, substitute(plain(k)[c] ~ ~vv, list(vv = format(values[nn], 
                ...))))
        }
        if (nn == "kg") {
            ans <- c(ans, substitute(plain(k)[g] ~ ~vv, list(vv = format(values[nn], 
                ...))))
        }
        if (nn == "s") {
            ans <- c(ans, substitute(s ~ ~vv, list(vv = format((values[nn]), 
                ...))))
        }
        if (nn == "logvarM") {
            ans <- c(ans, substitute({
                sigma^2
            }[M] ~ ~vv, list(vv = format(exp(values[nn]), ...))))
        }
        if (nn == "logNe") {
            ans <- c(ans, substitute(N[e] ~ ~vv, list(vv = format(exp(values[nn]), 
                ...))))
        }
        if (nn == "logn") {
            ans <- c(ans, substitute(n[e] ~ ~vv, list(vv = format(exp(values[nn]), 
                ...))))
        }
        if (nn == "o") {
            ans <- c(ans, substitute(theta ~ ~vv, list(vv = format(values[nn], 
                ...))))
        }
        if (nn == "relativekA0") {
            ans <- c(ans, substitute(k[A[0]] * "(%" * {
                sigma^2
            }[A[0]] * ")" ~ ~vv, list(vv = format(values[nn] * 
                100, ...))))
        }
        if (nn == "kA1") {
            ans <- c(ans, substitute(k[A[1]] ~ ~vv, list(vv = format(values[nn], 
                ...))))
        }
        if (nn == "kA2") {
            ans <- c(ans, substitute(k[A[2]] ~ ~vv, list(vv = format(values[nn], 
                ...))))
        }
        if (nn == "kA3") {
            ans <- c(ans, substitute(k[A[3]] ~ ~vv, list(vv = format(values[nn], 
                ...))))
        }
        if (nn == "relativekE0") {
            ans <- c(ans, substitute(k[E[0]] * "(%" * {
                sigma^2
            }[E[0]] * ")" ~ ~vv, list(vv = format(values[nn] * 
                100, ...))))
        }
        if (nn == "kE0") {
            ans <- c(ans, substitute(k[E[0]] ~ ~vv, list(vv = format(values[nn], 
                ...))))
        }
        if (nn == "kE1") {
            ans <- c(ans, substitute(k[E[1]] ~ ~vv, list(vv = format(values[nn], 
                ...))))
        }
        if (nn == "kE2") {
            ans <- c(ans, substitute(k[E[2]] ~ ~vv, list(vv = format(values[nn], 
                ...))))
        }
        if (nn == "kE3") {
            ans <- c(ans, substitute(k[E[3]] ~ ~vv, list(vv = format(values[nn], 
                ...))))
        }
        if (nn == "logrelativekA0") {
            ans <- c(ans, substitute(k[A[0]] * "(%" * {
                sigma^2
            }[A[0]] * ")" ~ ~vv, list(vv = format(exp(values[nn]) * 
                100, ...))))
        }
        if (nn == "logkA1") {
            ans <- c(ans, substitute(k[A[1]] ~ ~vv, list(vv = format(exp(values[nn]), 
                ...))))
        }
        if (nn == "logkA2") {
            ans <- c(ans, substitute(k[A[2]] ~ ~vv, list(vv = format(exp(values[nn]), 
                ...))))
        }
        if (nn == "logkA3") {
            ans <- c(ans, substitute(k[A[3]] ~ ~vv, list(vv = format(exp(values[nn]), 
                ...))))
        }
        if (nn == "logrelativekE0") {
            ans <- c(ans, substitute(k[E[0]] * "(%" * {
                sigma^2
            }[E[0]] * ")" ~ ~vv, list(vv = format(exp(values[nn]) * 
                100, ...))))
        }
        if (nn == "logkE0") {
            ans <- c(ans, substitute(k[E[0]] ~ ~vv, list(vv = format(exp(values[nn]), 
                ...))))
        }
        if (nn == "logkE1") {
            ans <- c(ans, substitute(k[E[1]] ~ ~vv, list(vv = format(exp(values[nn]), 
                ...))))
        }
        if (nn == "logkE2") {
            ans <- c(ans, substitute(k[E[2]] ~ ~vv, list(vv = format(exp(values[nn]), 
                ...))))
        }
        if (nn == "logkE3") {
            ans <- c(ans, substitute(k[E[3]] ~ ~vv, list(vv = format(exp(values[nn]), 
                ...))))
        }
    }
    if (!is.null(AIC)) {
        ans <- c(ans, substitute(AIC ~ ~vv, list(vv = format(AIC, 
            ...))))
    }
    return(ans)
}
