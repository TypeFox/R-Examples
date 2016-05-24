repnormmixmodel.sel <- function (x, k = 2, ...) 
{
    aic <- NULL
    bic <- NULL
    caic <- NULL
    icl <- NULL
        AIC <- function(emout) {
            emout$loglik - (length(emout$mu) + length(emout$stdev) + 
                length(emout$lambda) - 1)
        }
        BIC <- function(emout) {
            emout$loglik - log(nrow(x)) * (length(emout$mu) + 
                length(emout$stdev) + length(emout$lambda) - 
                1)/2
        }
        CAIC <- function(emout) {
            emout$loglik - (log(nrow(x)) + 1) * (length(emout$mu) + 
                length(emout$stdev) + length(emout$lambda) - 
                1)/2
        }
        ICL <- function(emout) {
            BIC(emout) - sum(emout$lambda * log(emout$lambda))
        }
        for (i in 1:k) {
            if (i == 1) {
                avx <- as.vector(x)
                mu <- mean(avx)
                s <- sd(avx)
                loglik <- sum(dnorm(avx, mean=mu, sd=s, log=TRUE))
                emout <- list(mu=mu, stdev=s, lambda=1, loglik=loglik)
            }
            else emout <- repnormmixEM(x, k = i, ...)
            aic[i] <- AIC(emout)
            bic[i] <- BIC(emout)
            caic[i] <- CAIC(emout)
            icl[i] <- ICL(emout)
        }
    out = rbind(aic, bic, caic, icl)
    Winner = apply(out, 1, function(x) (1:length(x))[x == max(x)])
    rownames(out) = c("AIC", "BIC", "CAIC", "ICL")
    colnames(out) = 1:k
    cbind(out, Winner)
}
