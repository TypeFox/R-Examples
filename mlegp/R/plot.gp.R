`plot.gp` <-
function(x, type = 0, params = NULL, sds = 1, CI.at.point=FALSE, ...) {

	if (type %in% c(2,3)) {
		if (is.null(params)) params = 1:x$numDim
		else params = toParamIndexes(params, x$params)
	}

	if (type <= 2) {
		if (sds > 0 && min(x$cv[,2]) < 0) {
			cat("warning: variance of predictor < 0 for at least one obs, probably due to roundoff error; will not plot standard deviation lines\n")
			sds = 0
		}

		if (sds > 0) {
			m1 = min(x$cv[,1] - sds*sqrt(x$cv[,2]))
			m1 = min(m1, x$Z)
			m2 = max(x$cv[,1] + sds*sqrt(x$cv[,2]))
			m2 = max(m2, x$Z)
		}
		else {
			m1 = min(x$cv[,1], x$Z)
			m2 = max(x$cv[,1], x$Z)
		}
	}

	
	if (type == 0) par(mfrow = c(2,1))
	else {
		if (type == 2 || type == 3) {
			## do not overwrite device if there is only 1 parameter
			if (length(params) != 1) createWindow(length(params))
		}
	}		

	if (type <=1) {
		plot(x$Z, x$cv[,1], ylim = c(m1, m2), xlab = "observed", ylab = "predicted", ...)
		index = order(x$Z)
		lines(x$Z[index], x$Z[index], lty=1)
		if (sds > 0) {
			if (!CI.at.point) {
				lines(x$Z[index], x$cv[index,1] - sds*sqrt(x$cv[index,2]), col = "red")
				lines(x$Z[index], x$cv[index,1] + sds*sqrt(x$cv[index,2]), col = "red")
			}
			else {
				for (i in 1:length(x$Z)) {
					lines(rep(x$Z[i],2),c(x$cv[i,1] - sds*sqrt(x$cv[i,2]), 
						x$cv[i,1] + sds *sqrt(x$cv[i,2])), col = "red"  )
				}
			}			
		}
	}

	if (type == 0) {
		plot(x$Z, t(x$Z - x$cv[,1]) / sqrt(x$cv[,2]), xlab="observed", 
		ylab = "standardized residual", ...)
		abline(0,0)
		abline(3,0, lty=9)
       	 	abline(-3,0, lty=9)
	}
	

	if (type == 2 || type == 3) {
		for (i in params) {
			if (type == 2) {
				plot(x$X[,i], x$cv[,1], ylim = c(m1, m2), xlab = x$param[i], ylab="predicted", ...)
				index = order(x$X[,i])
				lines(x$X[index,i], x$Z[index], lty=1)
				if (sds > 0) {
					if (!CI.at.point) {
					  lines(x$X[index,i], x$cv[index,1] - sds*sqrt(x$cv[index,2]), col = "red")
					  lines(x$X[index,i], x$cv[index,1] + sds*sqrt(x$cv[index,2]), col = "red")
					}
					else {
					  for (j in 1:length(x$X[,i])) {
						lines(rep(x$X[j,i],2),c(x$cv[j,1] - sds*sqrt(x$cv[j,2]), 
							x$cv[j,1] + sds *sqrt(x$cv[j,2])), 
						col = "red"  )
					  }
					}
				}
				
			}
			if (type == 3) {
		  		plot(x$X[,i], t(x$Z - x$cv[,1]) / sqrt(x$cv[,2]), 
					xlab=x$param[i], ylab = "standardized residual", ...)
	     	   		abline(0,0)
		       		abline(3,0, lty=9)
       			 	abline(-3,0, lty=9)
			}
		}
	}

	if (type == 4) {
	  z = t(x$Z - x$cv[,1]) / sqrt(x$cv[,2])
	  par(mfrow = c(1,2))
  	  qqnorm(z)
  	  qqline(z)
	  hist(z)       
	}

}

