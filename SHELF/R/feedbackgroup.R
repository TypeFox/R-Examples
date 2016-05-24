feedbackgroup <-
function(fit, quantiles =  NA, values = NA, dist = "best", sfg = 3){
	
	
	n.experts <- nrow(fit$limits)
	expertnames <- paste("expert.", 1:n.experts, sep="")
	
	distributions <- data.frame(matrix(0, nrow = 1, ncol = n.experts))
	names(distributions) <- expertnames
	distribution.names <- c("Normal", "Student-t", "Gamma", "Log normal", "Log Student-t", "Beta")
	

	
	if(is.na(quantiles[1]) == T ){
		quantiles <- fit$probs[1,]		
	}
	expert.quantiles <- data.frame(matrix(0, nrow = length(quantiles), ncol = n.experts), row.names = quantiles)
	names(expert.quantiles) <- expertnames
	
	if(is.na(values[1]) == T ){
		values <- fit$vals[1,]
	}
	expert.probs <- data.frame(matrix(0, nrow = length(values), ncol = n.experts), row.names = values)
	names(expert.probs) <- expertnames
	
	for(i in 1:n.experts){
		d.index<- c(1:2)
		if(fit$limits[i,1]>-Inf){
			d.index <- c(1:5)
			if(fit$limits[i,2]<Inf){
			d.index <- c(1:6)
			}
		}
		
		if(dist == "best"){
			d.select <- which(fit$ssq[i,] == min(fit$ssq[i, d.index]))[1]
		}else{
			d.select <- which(dist == c("normal", "t", "gamma", "lognormal", "logt","beta"))
		}
		
		distributions[1,i] <- distribution.names[d.select]
		
		temp <- feedbacksingle(fit, quantiles, values, ex = i)
		expert.quantiles[, i] <- temp$fitted.quantiles[, d.select]
		expert.probs[, i] <- temp$fitted.probabilities[, d.select]
	}
	
	
	list(expert.quantiles = signif(expert.quantiles, sfg), expert.probs = signif(expert.probs, sfg), distributions = distributions)
}
