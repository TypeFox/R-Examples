feedbacksingle <-
function(fit, quantiles =  NA, values = NA, sf = 3, ex = 1){
	
	n.distributions <- 6
	distribution.names <- c("Normal", "Student-t", "Gamma", "Log normal", "Log Student-t", "Beta")
	d.index<-c(1:2)
	
	if(is.na(quantiles[1]) == F ){
		report.elicited.q <- F
	}else{
		quantiles <- fit$probs[ex,]		
		report.elicited.q <- T	
	}
	
	Mq <- matrix(0, length(quantiles), n.distributions) 		
	Mq[,1] <- qnorm(quantiles, fit$Normal[ex,1], fit$Normal[ex,2])
	Mq[,2] <- qt(quantiles, fit$Student.t[ex,3]) * fit$Student.t[ex,2] + fit$Student.t[ex,1] 
	if(fit$limits[ex,1] > - Inf){
		d.index<-c(1:5)
		Mq[,3] <- fit$limits[ex,1] + qgamma(quantiles, fit$Gamma[ex,1], fit$Gamma[ex,2])
		Mq[,4] <- fit$limits[ex,1] + qlnorm(quantiles, fit$Log.normal[ex,1], fit$Log.normal[ex,2])
		Mq[,5] <- fit$limits[ex,1] + exp( qt(quantiles, fit$Log.Student.t[ex,3]) * fit$Log.Student.t[ex,2] + fit$Log.Student.t[ex,1])
		if(fit$limits[ex,2] < Inf){
			d.index<-c(1:6)
			Mq[,6] <- fit$limits[ex,1] + (fit$limits[ex,2] - fit$limits[ex,1]) * qbeta(quantiles, fit$Beta[ex,1], fit$Beta[ex,2] )
		}
	}
		
	if(is.na(values[1]) == F ){
		values <- matrix(values, nrow = length(values), ncol = n.distributions)
		report.elicited.p <- F
		}else{
		values <- matrix(fit$vals[ex,], nrow = length(fit$vals[ex,]), ncol = n.distributions)
		report.elicited.p <- T
		}
	
			
	values[,2] <- (values[,2] - fit$Student.t[ex,1]) / fit$Student.t[ex,2]
		
	if(fit$limits[ex,1] > - Inf){
		values[,3:4] <- values[,3:4] - fit$limits[ex,1]
		values[,5] <- (log(values[,5] - fit$limits[ex,1]) - fit$Log.Student.t[ex,1]) / fit$Log.Student.t[ex,2]
	}
		
	if((fit$limits[ex,1] > - Inf) & (fit$limits[ex,2] < Inf)){
		values[,6] <- (values[,6] - fit$limits[ex,1]) / (fit$limits[ex,2] - fit$limits[ex,1])
	}
		
	Mp <- matrix(0, nrow(values), ncol(values))
		
	Mp[,1] <- pnorm(values[,1], fit$Normal[ex,1], fit$Normal[ex,2])
	Mp[,2] <- pt(values[,2], fit$Student.t[ex,3])
	if(fit$limits[ex,1] > - Inf){
		Mp[,3] <- pgamma(values[,3], fit$Gamma[ex,1], fit$Gamma[ex,2])
		Mp[,4] <- plnorm(values[,4], fit$Log.normal[ex,1], fit$Log.normal[ex,2])
		Mp[,5] <- pt(values[,5], fit$Log.Student.t[ex,3])
		if(fit$limits[ex,2] <  Inf){
			Mp[,6] <- pbeta(values[,6], fit$Beta[ex,1], fit$Beta[ex,2])
		}
	}	
		
	if(report.elicited.p == F){
		Mp <- data.frame(Mp, row.names = values[,1])
		names(Mp) <- distribution.names
    dp.index<-d.index }else{
		Mp <- data.frame(matrix(fit$probs[ex,], ncol=1), Mp, row.names = values[,1])
		names(Mp) <- c("elicited", distribution.names)
    dp.index<-c(1,d.index+1)
	}
	
	if(report.elicited.q == F){
		Mq <- data.frame(Mq, row.names = quantiles)
		names(Mq) <- c(distribution.names)
		dq.index<-d.index }else{
		Mq <- data.frame(fit$vals[ex,], Mq, row.names = quantiles)
		names(Mq) <- c("elicited", distribution.names)
		dq.index<-c(1,d.index+1)
	}
	


	list(fitted.quantiles = signif(Mq[, dq.index],3), fitted.probabilities = signif(Mp[, dp.index],3))
}
