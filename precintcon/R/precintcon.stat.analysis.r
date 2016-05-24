#' @export
precintcon.stat.analysis <- function(..., args = NA) {
	
	l <- list(...)
	
	if (length(l) > 0) {
		
		data <- data.frame()
		
		pars <- ifelse(is.na(args), as.character(match.call()[1:length(l)+1]), args)
		
		for (i in 1:length(l)) {
			
			if (is.element("precintcon.daily",   class(l[[i]])) ||
				is.element("precintcon.monthly", class(l[[i]]))) {
				
				if (is.element("precintcon.monthly", class(l[[i]]))) { 
					
					data.mean.m  <- mean((l[[i]])[,3], na.rm=T)
					data.sd.m    <- sd((l[[i]])[,3],   na.rm=T)
					data.var.m   <- var((l[[i]])[,3],  na.rm=T)
					
					total        <- sum((l[[i]])[,3], na.rm=T)
					
					data <- rbind(data,
							data.frame(
									dataset=paste(pars[[i]], sep=""),
									mean.daily     = "---",
									sd.daily       = "---",
									var.daily      = "---",
									mean.monthly   = data.mean.m,
									sd.monthly     = data.sd.m,
									var.monthly    = data.var.m,
									total          = total								
							)
					)
				
				} else if (is.element("precintcon.daily", class(l[[i]]))) {
							
					data.mean    <- mean(as.vector(as.matrix((l[[i]])[,3:33])), na.rm=T)
					data.sd      <- sd(as.vector(as.matrix((l[[i]])[,3:33])), na.rm=T)
					data.var     <- var(as.vector(as.matrix((l[[i]])[,3:33])), na.rm=T)
						
					total        <- sum(as.vector(as.matrix((l[[i]])[,3:33])), na.rm=T)
					
					data.monthly <- precintcon.monthly.aggregation(l[[i]])
					
					data.mean.m  <- mean(data.monthly[,3], na.rm=T)
					data.sd.m    <- sd(data.monthly[,3],   na.rm=T)
					data.var.m   <- var(data.monthly[,3],  na.rm=T)
					
					data <- rbind(data, 
						data.frame(
							dataset=paste(pars[[i]], sep=""),
							mean.daily   = data.mean,
							sd.daily     = data.sd,
							var.daily    = data.var,
							mean.monthly = data.mean.m,
							sd.monthly   = data.sd.m,
							var.monthly  = data.var.m,
							total        = total
						))
							
				} 
			} else 
				stop("All input data should be either of classes \"precintcon.ci\", \"precintcon.fd\", \"precintcon.daily\", or \"precintcon.monthly\".")
		}
		
		ppp <- t(data)
		
		colnames(ppp) <- rep(" ", ncol(ppp))
				
		return(data)
	}
}
