print.topmod <-
function (x, ...) 
{
    tm = x
    if (length(tm$lik()) == 1) {
        infomat = c(tm$bool(), tm$lik(), tm$ncount())
        names(infomat) = c("Model Index", "Marg.Log.Lik.", "Sampled Freq.")
        print(infomat)
        betamat = cbind(as.vector(tm$betas_raw()), sqrt(as.vector(tm$betas2_raw()) - 
            as.vector(tm$betas_raw())^2))
        if (nrow(betamat) != 0) {
            if (ncol(betamat) == 1) {
                colnames(betamat) = "Coef."
            }
            else {
                colnames(betamat) = c("Coef.", "Std.Dev.")
            }
            rownames(betamat) = which(as.logical(as.vector(tm$bool_binary())))
            cat("\nEstimates:\n")
            print(betamat)
        }
        bin = as.vector(tm$bool_binary())
        names(bin) = 1:length(bin)
        cat("\nIncluded Covariates:\n")
        print(bin)
        cat("\nAdditional Statistics:\n")
        print(as.vector(tm$fixed_vector()))
    }
    else {
        mout = cbind(tm$lik(), tm$ncount())
        colnames(mout) = c("Marg.Log.Lik", "MCMC Freq")
        rownames(mout) = tm$bool()
        print(mout, ...)
    }
}
