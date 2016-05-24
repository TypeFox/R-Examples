ParamSampleCens <- function(censdata, dist = c("normal", "logistic", "gamma", "weibull")[1], null.values = c(0, 1), conf.level = 0.95, initial = NULL){

     censdata[is.na(censdata[, 1]) == FALSE & is.na(censdata[,2]) == TRUE, 2] <- Inf
     Delta <- matrix(1, nrow = nrow(censdata), ncol = 1)
     Delta[is.na(censdata[, 1]) == TRUE & is.na(censdata[, 2]) == FALSE] <- 0
     Delta[is.na(censdata[, 1]) == FALSE & is.na(censdata[, 2]) == FALSE & censdata[, 1] != censdata[, 2], 1] <- 0
     censdata[is.na(censdata[, 1]) == FALSE & is.na(censdata[, 2]) == FALSE & censdata[, 1] == censdata[, 2], 1] <- NA
     data <- censdata[, 2]
     LBInterval <- censdata[, 1]

     signif.level <- 1 - conf.level

     if(dist == "normal"){
        RowNames <- c("Mu", "Sigma")
        if (is.null(initial)){
               x <- apply(censdata, 1, mean, na.rm = TRUE)
               initial <- c(mean(x[x!=Inf]), sd(x[x!=Inf]))
        }
        MLE <- optim(par = initial, fn = LoglikNormalCens, data = data, lowerbound = LBInterval, vdelta = Delta, hessian = FALSE, control = list(fnscale = -1, reltol = 10^-15))
        HessianMatrix <- numDeriv::hessian(func = LoglikNormalCens, x = MLE$par, data = data, lowerbound = LBInterval, vdelta = Delta)
     }
    
     if(dist == "logistic"){
        RowNames <- c("Location", "Scale")
        if (is.null(initial)){initial <- c(0, 1)}
        MLE <- optim(par = initial, fn = LoglikLogisticCens, data = data, lowerbound = LBInterval, vdelta = Delta, hessian = FALSE, control = list(fnscale = -1))
        HessianMatrix <- numDeriv::hessian(func = LoglikLogisticCens, x = MLE$par, data = data, lowerbound = LBInterval, vdelta = Delta)
     }
    
     if(dist == "gamma"){
        RowNames <- c("Scale", "Rate")
        if (is.null(initial)){initial <- c(1, 1)}
        MLE <- optim(par = initial, fn = LoglikGammaCens, data = data, lowerbound = LBInterval, vdelta = Delta, hessian = FALSE, control = list(fnscale = -1))
        HessianMatrix <- numDeriv::hessian(func = LoglikGammaCens, x = MLE$par, data = data, lowerbound = LBInterval, vdelta = Delta)
     }
    
     if(dist == "weibull"){
        RowNames <- c("Shape", "Scale")
        if (is.null(initial)){initial <- c(1, 1)}
        MLE <- optim(par = initial, fn = LoglikWeibullCens, data = data, lowerbound = LBInterval, vdelta = Delta, hessian = FALSE, control = list(fnscale = -1))
        HessianMatrix <- numDeriv::hessian(func = LoglikWeibullCens, x = MLE$par, data = data, lowerbound = LBInterval, vdelta = Delta)
     }

     Estimation_param1 <- MLE$par[1]
     Estimation_param2 <- MLE$par[2]
	StandardError <- sqrt(-diag(solve(HessianMatrix)))
	CIparam1 <-  cbind(Estimation_param1 - qnorm(1 - signif.level / 2) * StandardError[1], Estimation_param1 + qnorm(1 - signif.level / 2) * StandardError[1])
	CIparam2 <-  cbind(Estimation_param2 - qnorm(1 - signif.level / 2) * StandardError[2], Estimation_param2 + qnorm(1 - signif.level / 2) * StandardError[2])

	percentage <- 1 - sum(Delta)/length(Delta)
	value_loglik <- MLE$value

     Results <- as.data.frame(matrix(nrow = 2, ncol = 5))
     colnames(Results) <- c("Estimator", "Std. Error", "CI.low", "CI.up", "p-value")
     rownames(Results) <- RowNames
     Results[1, 1] <- Estimation_param1
     Results[, 2] <- StandardError
	  Results[1, 3:4] <- CIparam1
	  Results[2, 1] <- Estimation_param2
	  Results[2, 3:4] <- CIparam2
     
     Results[, 1:4] <- round(Results[, 1:4], 5)

	p.value.test.param1 <- (Estimation_param1 - null.values[1]) / StandardError[1]
	p.value.param1 <- 1 - pchisq(p.value.test.param1 ^ 2, df = 1)
	p.value.test.param2 <- (Estimation_param2 - null.values[2]) / StandardError[2]
	p.value.param2 <- 1 - pchisq(p.value.test.param2 ^ 2, df = 1)
     Results[1, 5] <- format.pval(p.value.param1)
     Results[2, 5] <- format.pval(p.value.param2)
	     
	ReturnList <- list(coeff = Results,  percent.cens = percentage, loglik = value_loglik, info.converg = MLE$convergence, info.converg.message = MLE$message)

     return(ReturnList)
}








#
