`loglikfit` <-
function(angles, distribution, distpars=NA, ...){

		# Avoid zero densities:
		if(distribution %in% c("erectophile","plagiophile","extremophile")){
			angles[angles == 90] <- 90 - 1e-06
		}
		if(distribution == "extremophile"){
			angles[angles == 45] <- 45 - 1e-06
		}
		loglik <- sum(log(ftheta(angles, degrees=TRUE, distribution=distribution, distpars=distpars, ...)))
		
return(loglik)
}

