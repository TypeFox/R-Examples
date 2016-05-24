require("dse")

 data("eg1.DSE.data", package = "dse") 
 data("eg1.DSE.data.diff", package = "dse") 

  cat("truncate sample to 240 periods.\n")
  eg1.DSE.data.diff.trunc <- TSdata(input= inputData(eg1.DSE.data.diff)[1:240,, drop=F], 
            output=outputData(eg1.DSE.data.diff)[1:240,])
  seriesNames(eg1.DSE.data.diff.trunc)  <- seriesNames(eg1.DSE.data.diff)

  cat("estimates a VAR model using the truncated sample.\n")
  V.1 <- estVARXar(eg1.DSE.data.diff.trunc)
  
  cat("calculate the likelihood, one step ahead predictions, etc.\n")
  l.V.1 <-l(V.1, eg1.DSE.data.diff.trunc)
  
  cat("Likelihood and components for VAR model\n")
  # (with a breakdown for the 3 terms of the likelihood function).
  print(l.V.1$estimates$like, digits=16)

  #cat("Likelihood and components for VAR model\n")
  #l.V.1$estimates$like # also prints the value but not as many digits.

  cat("likelihood, one step ahead predictions, etc., based on the full sample.\n")
  o.V.1 <-l(V.1, eg1.DSE.data.diff) 

  cat("convert the VAR model to a state space model balanced by Mittnik's technique.\n")
  SS.V.1 <- toSS(V.1) 

  cat("likelihood, one step ahead predictions, etc., based on truncated sample.\n")
  l.SS.V.1 <-l(SS.V.1, eg1.DSE.data.diff.trunc) 

  cat("Likelihood and components for state space model\n")
  print(l.SS.V.1$estimates$like,digits=16) 

  cat("Maximum difference in one-step predictions of VAR and state space model ")
  # calculate the difference of the absolute values of the predictions of 
  #      the two models.
  cat(max(abs(l.V.1$estimates$pred - l.SS.V.1$estimates$pred))) 
  cat("\n")

  cat("Exhibit 2. Mittnik reduction from VAR model: \n")
  M5.SS.V.1 <- MittnikReduction(SS.V.1, data=eg1.DSE.data.diff, criterion="taic")  
  cat(paste(
   "  If criterion is not specified the program prompts for a state dimension\n",
   "  and returns that model. Results is put in the variable M5.SS.V.1."))

  cat("Exhibit 3. Mittnik estimation lag=3: \n")
  M12.shift3 <- estSSMittnik(eg1.DSE.data.diff.trunc, max.lag=3, n=12)
  M12.shift3 <- MittnikReduction(M12.shift3, data=eg1.DSE.data.diff.trunc, criterion="taic")  

  cat("Exhibit 4. Mittnik estimation lag=4: \n")
  M12.shift4 <- estSSMittnik(eg1.DSE.data.diff.trunc,max.lag=4, n=15)
  M12.shift4 <- MittnikReduction(M12.shift4, data=eg1.DSE.data.diff.trunc, criterion="taic")  


  cat(paste(
    "Plot cpi in year over year % change.\n",
    "Prediction is relative to previous month's actual (eg1.DSE.data)\n",
    "and % change is relative to actual.\n",
    "240 is the starting point for plotting.\n",
    "base is the start value of the undif, un logged series.\n"))
    
        i <- 3 # cpi is the third variable
	base <- eg1.DSE.data$output[1, i]
	pred <- o.V.1$estimates$pred[, i]
	y <- o.V.1$data$output[, i]
	y <- cumsum(c(log(base), y))
	pred <- c(log(base), pred)	# cumsum using pred relative to actual
	pred[2:length(pred)] <- pred[2:length(pred)] + y[1:(length(pred) - 1)]
	pred <- exp(pred)
	y <- exp(y)
	pred <- 100 * ((pred[13:length(pred)] - y[1:(length(y) - 12)])/y[1:(
		length(y) - 12)])
	y <- 100 * ((y[13:length(y)] - y[1:(length(y) - 12)])/y[1:(length(y) - 
		12)])
	tfplot(tfwindow(y, start=240),tfwindow(pred, start=240)) 

