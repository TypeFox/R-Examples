.dSigC <-
function(ndvi){
	days <- length(ndvi)
	ndvi <- .linIP(ndvi)
		
	doublesigmoid <- function(res) {
		##########single sigmoid#######
		sigmoid <- function(res, maxvalue=0, turnaround=FALSE) {
			days <- length(res)
			model <- vector(mode="numeric", length=days)
	
			if (turnaround){
				res <- res[days:1]
			}
		
			#calculate coefficients F and G
			#F - Base, (G+F) Maximum
			F <- mean(na.omit(res[1:(days/4)]))
			#F <- 0.2
			if (maxvalue > 0){
				G <- maxvalue - F
			} else {
				res.order <- order(na.omit(res), decreasing=TRUE)
				meanmax <- mean(res[res.order[1:3]])
				G <- meanmax - F
			}

			model <- .C("sigmoid", rdays=as.integer(days), ndvi=as.numeric(res),
				rF=as.numeric(F), rG=as.numeric(G), model=as.numeric(model), 
				PACKAGE="phenex")$model

			if (turnaround){
				model <- model[length(model):1]
			}
			if (length(model) != days){
				model <- rep(NA, days)
			}
			return(model)
		}

		#########################

		days <- length(res)
		delay <- 10
		count <- 1
		while((maxpos <- order(res, decreasing=TRUE)[count]) > 300){
			count <- count + 1
			if ((count > 100) || (is.na(maxpos))){
				return(rep(NA, days))
			}
		}
		#calculate front sigmoid
		modelfront <- sigmoid(res[1:(maxpos + delay)])
		if (length(na.omit(modelfront)) < 5){
			return(rep(NA, days))
		}
		#calculate back sigmoid
		maxvalue <- max(na.omit(modelfront))
		modelback <- sigmoid(res[(maxpos-1 + delay):days], maxvalue=maxvalue, turnaround=TRUE)
		if ((length(modelback) < 5) || is.na(modelback[3])){
			modelback <- rep(maxvalue, days-length(modelfront)+2)
		}	

		model <- vector(mode="numeric", length=days)

		#combine functions
		endfront <- length(modelfront)
		endback <- length(modelback)
		if (modelfront[endfront] != modelback[3]){
			scal <- modelback[3]
			for (i in 1:endback){
				modelback[i] <- (modelback[i] / scal) * modelfront[endfront]
			}
		}
		model <- c(modelfront, modelback[3:endback])
		
		return(model)
	}

	model <- doublesigmoid(ndvi)

	return(model)
}
