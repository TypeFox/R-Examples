require("dse")

 data("eg1.DSE.data", package = "dse") 
 data("eg1.DSE.data.diff", package = "dse") 

 cat("This demo reproduces some results from Bank of Canada Working Paper 93-4.\n")

   sub.sample <- TSdata(
      input=tfwindow(inputData(eg1.DSE.data.diff),end=c(1981,2)),
      output=tfwindow(outputData(eg1.DSE.data.diff),end=c(1981,2)) )

   VAR.model <- estVARXar(sub.sample, re.add.means=F)
   SS1.model <- l(balanceMittnik(toSS(VAR.model), n=9),sub.sample)

   g1 <- diag(1,9)
   g1[1:3,] <- SS1.model$model$H
   g1 <- solve(g1)
   g2 <- diag(1,9)
   g2[3,2:3] <- c(.1,-.1)
   g2[9,9] <- -1  # this is not really necessary but seems to have
                  #   happened in the paper
   example.gap.matrix <-g1 %*% g2
   SSgap.model <- l(gmap(example.gap.matrix,SS1.model),sub.sample)
   ARMA.model<- l(toARMA(SS1.model),sub.sample)

   print(VAR.model)   
   print(SS1.model)

   cat("Likelihood of VAR model:                          ")
   print(VAR.model$estimates$like[1], digits=16)
   cat("Likelihood of Mittnik balanced state space model: ")
   print(SS1.model$estimates$like[1], digits=16)
   cat("Likelihood of state space `gap' model:            ")
   print(SSgap.model$estimates$like[1], digits=16)
   cat("Likelihood of ARMA model:                         ")
   print(ARMA.model$estimates$like[1], digits=16)
   cat(paste(
    "Remark: A small change has been made in the likelihood\n",
    "calculation since the version of the code used for\n",
    "calculating the results in Bank of Canada Working Paper 93-4.\n",
    "The new method is more robust to degenerate densities but gives\n",
    "a small difference in the likelihood value. (The value reported \n",
    " was -2567.32801321424. )     P.Gilbert.\n"))
   
   cat("Stability of VAR model:\n")
   stability(VAR.model)
   
   cat("Stability of Mittnik balanced state space model:\n")
   stability(SS1.model)
   
   cat("Stability of state space `gap' model:\n")
   stability(SSgap.model)

   cat("Stability of ARMA model:\n")
   stability(ARMA.model)

   tfplot(VAR.model, Title="VAR model")
   
   cat(paste(
       "Remark: These are not advertised as best estimates. There is a bias.\n",
       "This estimation technique may be improved by setting some of the\n", 
       "options and other available estimation techniques work better.\n",
       "The example is intended primarily for illustrating the equivalence.\n"))

   tfplot(SS1.model, Title="Mittnik balanced state space model")

   tfplot(SSgap.model, Title="State space `gap' model")
   
   tfplot(ARMA.model,  Title="ARMA model")

   model<- l(VAR.model,eg1.DSE.data.diff)  # full sample


  cat(paste(
    "Plot cpi in year over year % change.\n",
    "Prediction is relative to previous month's actual (eg1.DSE.data)\n",
    "and % change is relative to actual.\n",
    "base is the start value of the undif, un logged series.\n"))

        i <- 3 # cpi is the third variable
	base <- eg1.DSE.data$output[1, i]
	pred <- model$estimates$pred[, i]
	y <- model$data$output[, i]
	y <- cumsum(c(log(base), y))
	pred <- c(log(base), pred)	# cumsum using pred relative to actual
	pred[2:length(pred)] <- pred[2:length(pred)] + y[1:(length(pred) - 1)]
	pred <- exp(pred)
	y <- exp(y)
	pred <- 100 * ((pred[13:length(pred)] - y[1:(length(y) - 12)])/y[1:(
		length(y) - 12)])
	y <- 100 * ((y[13:length(y)] - y[1:(length(y) - 12)])/y[1:(length(y) - 
		12)])

   cat("1 is the starting point for plotting.\n")
     tfplot(tfwindow(y, start=1),tfwindow(pred, start=1)) 
     title(main="Predicted and actual CPI in terms of per cent change \nover 12 months")

   cat("240 is the starting point for plotting.\n")
     tfplot(tfwindow(y, start=240),tfwindow(pred, start=240)) 
     title(main="Predicted and actual CPI in terms of per cent change \nover 12 months - ex post period")
