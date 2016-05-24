optspan=function (rdata, m = "adjusted", family = binomial(), verbose = TRUE, ...)
{
    if (tolower(m) == "adjusted") 
        x = length(rdata)
    else if (tolower(m) == "crude" | tolower(m) == "unadjusted") 
        x = 3
    else stop(paste("model type", m, "not recognized"))
    span.aic = matrix(ncol = 1, nrow = 19)
    for (S in 1:19) {
        if (x == 3) {
            fmla = as.formula(paste(names(rdata)[1], paste("lo(", 
                paste(names(rdata)[2:3], collapse = ","), ",span=", 
                S * 0.05, ")"), sep = "~"))
            fit = gam(fmla, family = family, data = rdata, ...)
        }
        else {
            fmla = as.formula(paste(names(rdata)[1], paste("lo(", 
                paste(names(rdata)[2:3], collapse = ","), ",span=", 
                S * 0.05, ")+", paste(names(rdata)[-(1:3)], collapse = "+")), 
                sep = "~"))
            fit = gam(fmla, family = family, data = rdata, ...)
        }
        span.aic[S, 1] = fit$aic
        if (verbose) cat(paste("The AIC for span=", format(S * 0.05, nsmall = 2), 
            " is ", format(fit$aic, nsmall = 2), ".", sep = ""), 
            fill = TRUE)
    }
    sp = which(span.aic[, 1] == min(span.aic[, 1])) * 0.05
    return(sp)
}
