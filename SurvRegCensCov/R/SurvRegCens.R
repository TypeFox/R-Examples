SurvRegCens <- function(formula, data = parent.frame(), Density, initial, conf.level = 0.95, 
                         intlimit = 10^-10, namCens = "VarCens", trace = 0, 
                         reltol = 10^-8){
     
     
     # ----------------------------------------------------------------------------
     # translate input arguments formula, data in the intially used input arguments
     #
     # initial version of function call was:
     #
     # SurvRegCens <- function(time, event, CovariateNonCens = NULL, 
     #                        CovariateCens, Density, initial, conf.level = 0.95, 
     #                        intlimit = 10^-10, namCens = "VarCens", trace = 0, 
     #                        reltol = 10^-8){     
     
     Call <- match.call()
     
     # pre-process input arguments
     indx <- match(c("formula", "data", "weights", "subset", "na.action"), names(Call), nomatch = 0)
     if (indx[1] == 0){stop("A formula argument is required")}
     temp <- Call[c(1, indx)]
     temp[[1]] <- as.name("model.frame") 
     m <- eval(temp, parent.frame())
     Terms <- attr(m, "terms")
     
     # extract and process response variable
     Y <- model.extract(m, component = "response")
     if (!inherits(Y, "Surv")){stop("Response must be a survival object")}
     time <- matrix(Y[, 1], ncol = 1)
     event <- as.vector(Y[, 2])
     
     # extract and process explanatory variables
     X <- model.frame(m)
     
     # identify where the interval-censored covariate is on RHS
     i.surv <- grep("Surv", colnames(X))
     
     # survival object on RHS
     r.surv <- as.matrix(X[, i.surv[-1]])
     CovariateCens <- matrix(NA, nrow = nrow(data), ncol = 2)
     
     # right-censored
     ind0 <- r.surv[, "status"] == 0
     CovariateCens[ind0, 1] <- r.surv[ind0, "time1"] 
     
     # events
     ind1 <- r.surv[, "status"] == 1
     CovariateCens[ind1, 1:2] <- r.surv[ind1, "time1"] 
     
     # left-censored
     ind2 <- r.surv[, "status"] == 2
     CovariateCens[ind2, 2] <- r.surv[ind2, "time1"] 
     
     # left-censored
     ind3 <- r.surv[, "status"] == 3
     CovariateCens[ind3, 1] <- r.surv[ind3, "time1"] 
     CovariateCens[ind3, 2] <- r.surv[ind3, "time2"] 
     
     # extract further explanatory variables used in formula on RHS, if any
     nam1 <- colnames(X)[-i.surv]
     CovariateNonCens <- X[, -i.surv]
     if (length(colnames(X)[-i.surv]) == 0){CovariateNonCens <- NULL}
     if (is.null(CovariateNonCens) == FALSE){CovariateNonCens <- as.matrix(CovariateNonCens)}
     # ----------------------------------------------------------------------------
     
     signif.level <- 1 - conf.level
     
     CovariateCens[is.na(CovariateCens[, 1]) == FALSE & is.na(CovariateCens[,2]) == TRUE, 2] <- Inf
     CovariateCens2 <- CovariateCens
     CovariateCens2[is.na(CovariateCens[, 1]) == FALSE & is.na(CovariateCens[, 2]) == FALSE & CovariateCens[, 1] == CovariateCens[, 2], 1] <- NA
     VectorR <- matrix(1, nrow = nrow(CovariateCens), ncol = 1)
     VectorR[is.na(CovariateCens[, 1]) == TRUE & is.na(CovariateCens[, 2]) == FALSE] <- 0
     VectorR[is.na(CovariateCens[, 1]) == FALSE & is.na(CovariateCens[, 2]) == FALSE & CovariateCens[, 1] != CovariateCens[, 2]] <- 0
     VectorLB <- CovariateCens2[, 1]
     
     # do maximization
     result_value <- optim(initial, LoglikWeibullSurvRegCens, data_y = time, data_cov_noncens = CovariateNonCens, 
                           data_cov_cens = CovariateCens[, 2], density = Density, data_r_loglik = VectorR, 
                           data_lowerbound = VectorLB, data_delta_loglik = event, intlimit = intlimit, 
                           control = list(maxit = 5000, fnscale = -1, trace = trace, reltol = reltol), hessian = FALSE)
     
     if(is.null(CovariateNonCens) == FALSE){
          NumberNonCens <- ncol(as.matrix(CovariateNonCens))
          NamesBeta <- matrix(nrow = 1, ncol = NumberNonCens)
     }
     
     Estimation_lambda <- result_value$par[1]
     Estimation_gamma <- result_value$par[2]
     if(is.null(CovariateNonCens) == FALSE){Estimation_betaNonCens <- result_value$par[3:(length(result_value$par) - 1)]}
     Estimation_betaCens <- result_value$par[length(result_value$par)]
     
     # standard errors
     HessianMatrix <- numDeriv::hessian(func = LoglikWeibullSurvRegCens, x = result_value$par, data_y = time, 
                                        data_cov_noncens = CovariateNonCens, data_cov_cens = CovariateCens[, 2], 
                                        density = Density, data_r_loglik = VectorR, data_lowerbound = VectorLB, 
                                        data_delta_loglik = event, intlimit = intlimit)
     SEs <- sqrt(-diag(solve(HessianMatrix)))
     
     qa <- qnorm(1 - signif.level / 2)
     CIlambda <-  cbind(Estimation_lambda - qa * SEs[1], Estimation_lambda + qa * SEs[1])
     CIgamma <-  cbind(Estimation_gamma - qa * SEs[2], Estimation_gamma + qa * SEs[2])
     
     if(is.null(CovariateNonCens) == FALSE){
          for(ii in 1:NumberNonCens){
               assign(paste("CIBetaNonCens", ii, sep = ""), cbind(Estimation_betaNonCens[ii] - qa * SEs[ii+2], 
                                                                  Estimation_betaNonCens[ii] + qa * SEs[ii+2]))
               NamesBeta[ii] <- paste("Beta", ii, sep = "")
          }
     }
     
     CIbetaCens <- cbind(Estimation_betaCens - qa * SEs[length(SEs)], Estimation_betaCens + qa * SEs[length(SEs)])
     
     if(is.null(CovariateNonCens) == FALSE){
          CIbetaNonCens <- matrix(nrow = 2 ,ncol = NumberNonCens)
          for(ii in 1:NumberNonCens){
               CIbetaNonCens[,ii] <- get(paste("CIBetaNonCens", ii, sep = ""))
          }
          
          if(NumberNonCens >= 2){
               Estimation_betaNonCens <- t(as.matrix(Estimation_betaNonCens))
               colnames(Estimation_betaNonCens) <- NamesBeta
               colnames(CIbetaNonCens) <- NamesBeta
          }
     }
     
     # collect results
     table <- as.data.frame(matrix(nrow = length(result_value$par), ncol = 8))
     tNamesRow <- matrix(nrow = nrow(table), ncol = 1)
     tNamesRow[1] <- "lambda"
     tNamesRow[2] <- "gamma"
     if(is.null(CovariateNonCens) == FALSE){
          if(ncol(as.matrix(CovariateNonCens)) > 1){
               tNamesRow[3:(nrow(tNamesRow) - 1)] <- nam1  ## paste("BetaNonCens", seq(from = 1, to = ncol(as.matrix(CovariateNonCens)), by = 1), sep = "")
          }          
          if(ncol(as.matrix(CovariateNonCens)) == 1){tNamesRow[3] <- nam1}
     } 
     
     tNamesRow[nrow(tNamesRow)] <- namCens
     rownames(table) <- tNamesRow
     colnames(table) <- c("Estimate", "Std. Error", "CI.low", "CI.up", "p-value", "exp(Estimator)", "exp(CI.low)", "exp(CI.up)")
     
     table[, 1] <- result_value$par
     table[, 2] <- SEs
     table[1, 3:4] <- CIlambda
     table[2, 3:4] <- CIgamma
     
     if(is.null(CovariateNonCens) == FALSE){table[3:(2 + ncol(CIbetaNonCens)), 3:4] <- t(CIbetaNonCens)}
     table[nrow(table), 3:4] <- CIbetaCens
     table[3:nrow(table), 6] <- exp(table[3:nrow(table), 1])
     table[3:nrow(table), 7:8] <- exp(table[3:nrow(table), 3:4])
     
     if(is.null(CovariateNonCens) == FALSE){
          p.value.BetaNonCens <- matrix(nrow = 1, ncol = NumberNonCens)
          for(ii in 1:NumberNonCens){
               test.BetaNonCens_ii <- (Estimation_betaNonCens[ii] - 0) / SEs[ii + 2]
               p.value.BetaNonCens[ii] <- 1 - pchisq(test.BetaNonCens_ii ^ 2, df = 1)
          }                                                                
          
          if(NumberNonCens >= 2){colnames(p.value.BetaNonCens) <- NamesBeta}
     }
     
     test.BetaCens <- (Estimation_betaCens - 0) / SEs[length(SEs)]
     p.value.BetaCens <- 1 - pchisq(test.BetaCens ^ 2, df = 1)
     if(is.null(CovariateNonCens) == FALSE){table[3:(nrow(table) - 1), 5] <- as.vector(p.value.BetaNonCens)}
     
     table[nrow(table), 5] <- p.value.BetaCens     
     percentage <- (1 - (sum(VectorR) / length(VectorR))) * 100
     loglik <- result_value$value
     
     # AIC, BIC
     n <- length(time)
     d <- nrow(table)
     AIC <- - 2 * loglik + 2 * d
     BIC <- - 2 * loglik + d * log(n)
          
     res <- list("coeff" = table, "percent.cens" = percentage, "loglik" = loglik, "d" = d, "n" = n,
                        "AIC" = AIC, "BIC" = BIC, "info.converg" = result_value$convergence, 
                        "info.converg.message" = result_value$message, "Call" = Call)
     class(res) <- "src"
     return(res)
}





#