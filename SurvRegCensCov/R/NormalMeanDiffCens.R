NormalMeanDiffCens <- function(censdata1, censdata2, conf.level = 0.95, null.values = c(0, 0, 1, 1)){
    
     # initial=(mu1, delta, sigma1, sigma2); delta = mu1 - mu2
     initial <- c(mean(censdata1[is.na(censdata1[, 2]) == FALSE, 2]), mean(censdata1[is.na(censdata1[, 2]) == FALSE, 2]) - mean(censdata2[is.na(censdata2[, 2]) == FALSE, 2]), sd(censdata1[is.na(censdata1[, 2]) == FALSE, 2]), sd(censdata2[is.na(censdata2[, 2]) == FALSE, 2]))
     signif.level <- 1 - conf.level
     
     censdata1[is.na(censdata1[, 1]) == FALSE & is.na(censdata1[, 2]) == TRUE, 2] <- Inf
     Delta1 <- matrix(1, nrow = nrow(censdata1), ncol = 1)
     Delta1[is.na(censdata1[, 1]) == TRUE & is.na(censdata1[, 2]) == FALSE] <- 0
     Delta1[is.na(censdata1[, 1]) == FALSE & is.na(censdata1[, 2]) == FALSE & censdata1[, 1] != censdata1[, 2], 1] <- 0
     data1 <- censdata1[, 2]
     LBInterval1 <- censdata1[, 1]
     LBInterval1[is.na(censdata1[, 1]) == FALSE & is.na(censdata1[, 2]) == FALSE & censdata1[, 1] == censdata1[, 2]] <- NA
     
     censdata2[is.na(censdata2[,1])==FALSE & is.na(censdata2[,2])==TRUE,2] <- Inf
     Delta2 <- matrix(1, nrow = nrow(censdata2), ncol = 1)
     Delta2[is.na(censdata2[, 1]) == TRUE & is.na(censdata2[, 2]) == FALSE] <- 0
     Delta2[is.na(censdata2[, 1]) == FALSE & is.na(censdata2[, 2]) == FALSE & censdata2[, 1] != censdata2[, 2], 1] <- 0
     data2 <- censdata2[, 2]
     LBInterval2 <- censdata2[, 1]
     LBInterval2[is.na(censdata2[, 1]) == FALSE & is.na(censdata2[, 2]) == FALSE & censdata2[, 1] == censdata2[, 2]] <- NA

	    MLE <- optim(par = initial, fn = LoglikNormalDeltaCens, data1 = data1, lowerbound1 = LBInterval1, vdelta1=Delta1, data2 = data2, lowerbound2 = LBInterval2, vdelta2 = Delta2, hessian = FALSE, control = list(fnscale = -1, reltol = 10^-15))
  	  Estimation_delta <- MLE$par[2]
     
     
     ## compute hessian using numDeriv     
     HessianMatrix <- numDeriv::hessian(func = LoglikNormalDeltaCens, x = MLE$par, data1 = data1, lowerbound1 = LBInterval1, vdelta1=Delta1, data2 = data2, lowerbound2 = LBInterval2, vdelta2 = Delta2)
	   StandardError <- sqrt(-diag(solve(HessianMatrix)))
	    CIdelta <-  cbind(Estimation_delta - qnorm(1 - signif.level / 2) * StandardError[2], Estimation_delta + qnorm(1 - signif.level / 2) * StandardError[2])
    
	    p.value.test <- Estimation_delta / StandardError[2]
	    p.value <- 1 - pchisq(p.value.test ^ 2, df = 1)
    
     res1 <- data.frame(matrix(NA, ncol = 5, nrow = 5))
     colnames(res1) <- c("Estimator", "Std. Error", "CI.low", "CI.up", "p-value")
     rownames(res1) <- c("mu1", "mu2", "sigma1", "sigma2", "Mean difference delta")
     
     res1[5, 1] <- Estimation_delta
     res1[5, 2] <- StandardError[2]
     res1[5, 3:4] <- CIdelta
     res1[5, 5] <- format.pval(p.value)
    
     MLE1 <- ParamSampleCens(censdata = censdata1, null.values = null.values[c(1, 3)], initial = c(initial[1], initial[3]))$coeff
     res1[c(1, 3), 1:5] <- MLE1
     MLE2 <- ParamSampleCens(censdata = censdata2, null.values = null.values[c(2, 4)], initial = c(initial[1] - initial[2], initial[4]))$coeff
     res1[c(2, 4), 1:5] <- MLE2
         
     return(res1)
}









#
