# New version of this -- optimized for speed and more values of the
# hyperparameters

SZ.prior.evaluation <- function (Y, p, lambda0, lambda1, lambda3, lambda4, lambda5,
                                 mu5, mu6, z = NULL, nu = ncol(Y) + 1,
                                 qm,
                                 prior = 0, nsteps,
                                 y.future)
{
    combos <- length(lambda0) * length(lambda1) * length(lambda3) *
        length(lambda4) * length(lambda5) * length(mu5) * length(mu6)

    results <- matrix(0, combos, 11)

    results[,1:7] <- as.matrix(expand.grid(lambda0=lambda0, lambda1=lambda1,
                                           lambda3=lambda3, lambda4=lambda4,
                                           lambda5=lambda5, mu5=mu5, mu6=mu6))

    for (i in 1:nrow(results))
    {
        fit <- szbvar(Y, p, z, lambda0=results[i,1],
                      lambda1=results[i,2],
                      lambda3=results[i,3], lambda4=results[i,4],
                      lambda5 = results[i,5],
                      mu5 = results[i,6], mu6 = results[i,7], nu =
                      ncol(Y)+1, qm = qm, prior = prior,
                      posterior.fit = T)

        forecast <- forecast(fit, nsteps)
        eval.forecasts <- cf.forecasts(forecast[(nrow(forecast) -
                                                 nsteps + 1):nrow(forecast), ], y.future)

        tmp <- c(eval.forecasts[1], eval.forecasts[2], fit$marg.llf[1], fit$marg.post[1])

        results[i, 8:11] <- tmp

        if(i%%100==T) { cat("Finished model", i, "of", nrow(results),"\n") }

    }
    colnames(results) <- c("lambda0", "lambda1", "lambda3", "lambda4",
                           "lambda5", "mu5", "mu6", "RMSE", "MAE", "LLF", "logMDD")
    return(results)
}
