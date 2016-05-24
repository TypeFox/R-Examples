#' Print a summary of the fitted model.

#' @return The posterior mean and 95 percent credible intervals, n_eff, Rhat and WAIC.
#' @param x A cdtafit object from \link{fit}.
#' @param digits An optional positive value to control the number of digits to print when printing numeric values. The default is 3.
#' @param ... other \link[rstan]{stan} options.
#' @examples
#'
#' \dontrun{
#'
#' fit1 <- fit(data=telomerase,
#'              SID = "ID",
#'              copula="fgm",
#'              iter = 400,
#'              warmup = 100,
#'              seed=1,
#'              cores=1)
#'
#' print(fit1)
#'
#'}
#' @references {Watanabe S (2010). Asymptotic Equivalence of Bayes Cross Validation and Widely Applicable Information Criterion in Singular
#' Learning Theory. Journal of Machine Learning Research, 11, 3571-3594.}
#' @references {Vehtari A, Gelman A (2014). WAIC and Cross-validation in Stan. Unpublished, pp. 1-14.}
#' @export
#' @author Victoria N Nyaga \email{victoria.nyaga@outlook.com}
print.cdtafit <- function(x, digits=3, ...){

#=======================Extract Model Parameters ===================================#
   sm <- summary.cdtafit(x, ...)

   mu <- data.frame(sm$allsm$summary[grepl('MU', rownames(sm$allsm$summary)), c("mean", "2.5%", "97.5%", "n_eff", "Rhat")])

    if (nrow(mu) > 2){
        ktau <- data.frame(sm$allsm$summary[grepl('ktau', rownames(sm$allsm$summary)), c("mean", "2.5%", "97.5%", "n_eff", "Rhat")])
        rr <- data.frame(sm$allsm$summary[grepl('RR', rownames(sm$allsm$summary)), c("mean", "2.5%", "97.5%", "n_eff", "Rhat")])
    } else {
        ktau <- sm$allsm$summary[grepl('ktau', rownames(sm$allsm$summary)), c("mean", "2.5%", "97.5%", "n_eff", "Rhat")]
    }
#==========================Tranform omega to ktau in FRANK =========================================#
    if (x@copula=="frank"){
        if (nrow(mu) > 2){
            omega <- data.frame(sm$allsm$summary[grepl('betaomega', rownames(sm$allsm$summary)), c("mean", "2.5%", "97.5%", "n_eff", "Rhat")])
            for(i in 1:nrow(ktau)){
                ktau[i,1] <- omega.to.ktau(omega[i,1])
                ktau[i,2] <- omega.to.ktau(omega[i,2])
                ktau[i,3] <- omega.to.ktau(omega[i,3])
            }
        } else {
            omega <- sm$allsm$summary[grepl('betaomega', rownames(sm$allsm$summary)), c("mean", "2.5%", "97.5%", "n_eff", "Rhat")]
            ktau[1] <- omega.to.ktau(omega[1])
            ktau[2] <- omega.to.ktau(omega[2])
            ktau[3] <- omega.to.ktau(omega[3])
        }
    }
#===================================    =======         ============================================#
    if (nrow(mu) > 2){
        Summary <- rbind(mu, rr, ktau)
    } else {
        Summary <- rbind(mu, ktau)
        row.names(Summary)[3] <- "ktau[1]"
    }
#========================== ============================= =========================================#
    if (nrow(mu) > 2){
        Summary$Parameter <- c(rep(c("Sensitivity", "Specificity"), each=nrow(mu)/2),
                               rep(c("Sensitivity", "Specificity"), each=nrow(rr)/2),
                               rep("Correlation", each=nrow(ktau)))
    } else {

        Summary$Parameter <- c(rep(c("Sensitivity", "Specificity"), each=nrow(mu)/2), "Correlation")
    }

    names(Summary) <- c("Mean", "Lower", "Upper", "n_eff", "Rhat", "Parameter")

    Summary <- Summary[,c(6, 1:5)]
    cat("Posterior marginal mean sensitivity and specificity\n\twith 95% credible intervals\n\n")
    print(Summary, digits=digits)
    cat("\n\n")

    cat("Model characteristics\n\n")
    cat(paste("Copula function: ", x@copula,  sep=""))
    cat(paste(", sampling algorithm: ", attr(x@fit@sim$samples[[1]], "args")$sampler_t, "\n", sep=""))
    cat(paste("\nFormula(1):  MUse ~ ", as.character(x@modelargs$formula.se)[3], sep=""))
    cat(paste("\nFormula(2):  MUsp ~ ", as.character(x@modelargs$formula.sp)[3], sep=""))
    cat(paste("\nFormula(3):  Omega ~ ", as.character(x@modelargs$formula.omega)[3], sep=""))


    cat(paste("\n", x@fit@sim$chains, " chain(s)", "each with iter=", x@fit@sim$iter,"; ", "warm-up=",
              x@fit@sim$warmup, "; ", "thin=", x@fit@sim$thin, ".\n",  sep=""))
    cat(paste("post-warmup draws per chain=",
              (x@fit@sim$iter-x@fit@sim$warmup)/x@fit@sim$thin, ";", "total post-warmup draws=",
              ((x@fit@sim$iter-x@fit@sim$warmup)/x@fit@sim$thin)*x@fit@sim$chains, ".\n", sep=""))

    w <- waic(x@fit)

    cat("\nPredictive accuracy of the model\n\n")
    cat(paste("Log point-wise predictive density (LPPD): ", sprintf(paste("%.", digits, "f", sep=''),w$lppd), sep=''))
    cat("\n")
    cat(paste("Effective number of parameters: ",  sprintf(paste("%.", digits, "f", sep=''),w$p_waic), sep=''))
    cat("\n")
    cat(paste("Watanabe-Akaike information Criterion (WAIC): ", sprintf(paste("%.", digits, "f", sep=''),w$waic), sep=''))
    cat("\n\n")
}

