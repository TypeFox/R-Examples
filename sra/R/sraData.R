sraData <-
function (phen.mean, phen.var, phen.sel, var.sel = NULL, N = NULL, 
    gen = NULL, rep = NULL) 
{
    "vartruncGauss" <- function(mu, sigma, mustar) {
        if (any(c(is.na(mu), is.na(mustar)))) {
            return(sigma * sigma)
        }
        if (mu > mustar) {
            mustar <- mu + mu - mustar
        }
        if (isTRUE(all.equal(mu, mustar))) {
            return(sigma * sigma)
        }
        ff <- function(alpha) {
            dnorm(alpha)/(1 - pnorm(alpha))
        }
        lambda <- (mustar - mu)/sigma
        alpha <- uniroot(function(x) {
            ff(x) - lambda
        }, lower = -10, upper = 10)$root
        return(sigma * sigma * (1 - lambda * (lambda - alpha)))
    }
    if (is.null(gen)) {
        gen <- 1:(length(phen.mean))
    }
    if (is.null(rep)) {
        for (g in 1:length(gen)) {
            rep <- c(rep, ifelse(g > 1, 1 + sum(gen[1:(g - 1)] == 
                gen[g]), 1))
            if (max(rep) > 1) {
                warning("rep should be provided; ambiguous data set. Results are error-prone.")
            }
        }
    }
    if (is.null(N)) {
        N <- rep(100, length(phen.mean))
        warning("Population size not provided, considered as N=100.")
    }
    rep <- as.factor(rep)
    if (length(phen.mean) != length(phen.var)) {
        stop("Not the same number of means and variances. Use NA if necessary.")
    }
    if (length(phen.mean) == length(phen.sel) - 1) {
        phen.sel <- c(phen.sel, NA)
    }
    if (length(phen.mean) != length(phen.sel)) {
        stop("Not the same number of means and selected parents. Use NA if necessary.")
    }
    if (length(phen.mean) != length(gen)) {
        stop("The number of generations does not match the data. Use NA if necessary.")
    }
    if (length(phen.mean) != length(rep)) {
        stop("The number of repetitions does not match the data.")
    }
    if (!is.null(var.sel) && (length(phen.mean) != length(var.sel))) {
        stop("The number of generations does not match the data. Use NA if necessary.")
    }
    if (is.null(var.sel)) {
        warning("The phenotypic variance of breeders was not provided.\n An educated guess will be made, model fitting might be inaccurate.\n")
        for (i in 1:length(phen.mean)) {
            var.sel <- c(var.sel, vartruncGauss(mu = phen.mean[i], 
                sigma = sqrt(phen.var[i]), mustar = phen.sel[i]))
        }
    }
    ans <- data.frame(rep = rep, gen = gen, mean = phen.mean, 
        var = phen.var, sel = phen.sel, vsel = var.sel, N = N)
    class(ans) <- c("sradata", class(ans))
    return(ans)
}
