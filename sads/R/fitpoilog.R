fitpoilog <- function(x, trunc = 0, ...){
    dots <- list(...)
    if (!is.null(trunc)){
        if (min(x)<=trunc) stop("truncation point should be lower than the lowest data value")
        else{
            if(trunc==0){
                pl.par <- poilogMLE(x, startVals = c(mu = mean(log(x)) + log(0.5), sig = sd(log(x))), zTrunc = TRUE)$par
            }
            else pl.par <- poilogMLE(x, startVals = c(mu = mean(log(x)) + log(0.5), sig = sd(log(x))))$par
            LL <- function(mu, sig) -sum(dtrunc("poilog", x = x, coef = list(mu = mu, sig = sig), trunc = trunc, log = TRUE))
        }
    }
    if (is.null(trunc)){
        pl.par <- poilogMLE(x, startVals = c(mu = mean(log(x)) + log(0.5), sig = sd(log(x))), zTrunc = FALSE)$par
        LL <- function(mu, sig) -sum(dpoilog(x, mu, sig, log = TRUE))
    }
    result <- do.call("mle2", c(list(LL, start = as.list(pl.par), data = list(x = x)), dots))
    new("fitsad", result, sad="poilog", distr = distr.depr, trunc = ifelse(is.null(trunc), NaN, trunc))
}
