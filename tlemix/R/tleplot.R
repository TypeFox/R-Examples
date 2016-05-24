setGeneric("tleplot",
		function(object,data, ...)
			standardGeneric("tleplot"))


setMethod("tleplot",
		signature(object="TLE",data="data.frame"),	
		function(object, data,...) {
			nobs <- dim(data)[1] # number of observations
			outliers <- object@indout # index of outliers
			normalobs <- object@indbest # "normal" observations
			col <- rep(1,nobs) # create initial color vector
			
			
			# set color to clusternumber of the observations
			# that where selected by flexmix
			#col[normalobs] <- object@estimate@cluster+1
			col[normalobs] <- object@tleclusters[object@indbest]+1
			
			## assert: col is a vector of length nobs containing
			# ones (=black) for outliers, and > 1 for normal observations
			dottype <- rep(16,nobs)
			dottype[outliers] <- rep(2,length(outliers))
			# draw a scatterplot for twodimensional
			#### CAVE, TODO: x,y selected by hand, needs to be generalized
			plot(data$x,data$y,col=col,pch=dottype,...)
			
			## Legend
			cf <- as.factor(object@estimate@cluster)
			n <- length(levels(cf))
			string <- character(length=n+1)
			for(i in 1:n) {
				string[i] <- paste("Cluster", i)
			}
			string[n+1] <- "Outliers"
			legend("bottomright", string, col=(n+1):1, pch=cbind(rep(16,n),2))
		})
