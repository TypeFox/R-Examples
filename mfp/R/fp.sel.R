fp.sel <- function(fit, alpha = 0.05, select = 1)
{
#
# Calculate deviance differences & p-values
#
dispersion <- fit$dispersion 
#
# 4 d.f. test at the alpha level of the best-fitting FP2 against the null model
    dd.null <- fit$dev0 - min(fit$dev1, fit$dev2, fit$dev4, na.rm=TRUE)  
    p.null <- pchisq(dd.null/dispersion, fit$df, lower.tail=FALSE)
    if(fit$df > 1) {
# 3 d.f. test at the alpha level of the best-fitting FP2 against a straight line
        dd.lin <- fit$dev1 - min(fit$dev2, fit$dev4, na.rm=TRUE)
        p.lin <- pchisq(dd.lin/dispersion, fit$df - 1, lower.tail=FALSE)
        if(fit$df > 2) {
# 2 d.f. test at the alpha level of the best-fitting FP2 against the best-fitting FP1
            dd.FP <- fit$dev2 - fit$dev4
            p.FP <- pchisq(dd.FP/dispersion, 2, lower.tail=FALSE)
        }
        else p.FP <- NA
    }
    else {
        p.lin <- NA
        p.FP <- NA
    }
    if(p.null > select) {
        df <- 0
        pwrs <- c(NA, NA)
        dev <- fit$dev0
    }
    else {
        if(fit$df > 1) {
            if(p.lin > alpha)
                df <- 1
            else {
                if(fit$df > 2) {
                  if(p.FP > alpha)
                    df <- 2
                  else {
                    df <- 4
                    pwrs <- fit$pwr4
                    dev <- fit$dev4
                  }
                }
                else df <- 2
            }
        }
        else df <- 1
    }
    if(df == 1) {
        pwrs <- c(1, NA)
        dev <- fit$dev1
    }
    if(df == 2) {
        pwrs <- c(fit$pwr2, NA)
        dev <- fit$dev2
    }
    results <- list(p.null = p.null, p.lin = p.lin, p.FP = p.FP, df = df, 
        pwrs = pwrs, dev = dev)
    return(list(results=results, fit=fit))
}
