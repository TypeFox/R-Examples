 conf.int.lm <- function (dat.lm, alpha){ 

	# variables and names  
	 lm.x <- dat.lm$model[,2] 		# extract x variable
         lm.y <- dat.lm$model[,1] 		# extract y variable
         name.x <- names(dat.lm$model[2])	# name of x
         name.y <- names(dat.lm$model[1])	# name of y  

	# calculate stderr and conf int
         rmse <- sqrt(sum(dat.lm$resid^2)/dat.lm$df) # calculates residual standard error
         more <- ((lm.x - mean(lm.x))^2)/var(lm.x) # square dev over variance
         stderr <- (rmse/sqrt(dat.lm$df+2))*(1+ more) # std error of estimates 
         t.value <- qt(1 - alpha/2, dat.lm$df) # calculates t value for given alpha
         lower <- dat.lm$fitted - stderr*t.value # confidence interval low end
         upper <- dat.lm$fitted + stderr*t.value # confidence interval high end

	# graphics
         plot(lm.x, lm.y, xlab=name.x, ylab=name.y)	# scatter plot 
         abline(dat.lm$coef)	# regression line
         ord <- order(lm.x)	# sort 
         lines(lm.x[ord], lower[ord])	# plot lower
         lines(lm.x[ord], upper[ord])	# plot high
         identify(lm.x,lm.y,labels=row.names(dat.lm$model)) 
	# return
	 invisible(list(lower = lower, upper = upper))
 }
