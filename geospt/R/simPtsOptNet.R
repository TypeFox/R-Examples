simPtsOptNet <- function(formula, loc=NULL, data, fitmodel, BLUE=FALSE, n, popSize, generations, xmin, ymin, xmax, ymax, plotMap=FALSE, spMap=NULL, ...){
	evaluate <- function(string=c()) {
		returnVal = NA;
		pts2 <- as.data.frame(matrix(0, ncol=ncol(as.data.frame(data)), nrow=n))
                names(pts2) <- colnames(as.data.frame(data))
                
                if(is.data.frame(data)) {
                if(is.null(loc)) stop(paste("loc must be provided"))
                x1 <- all.vars(loc)[1]
                y1 <- all.vars(loc)[2]
		}
                
                if(class(data)=="SpatialPointsDataFrame") {
                x1 <- colnames(coordinates(data))[1]
                y1 <- colnames(coordinates(data))[2]
		}

                for (i in 1:n){
			pts2[i,x1] <- round(string[i], 1)
		}
		for (j in 1:n){
			pts2[j,y1] <- round(string[n + j], 1)
		}
                coordinates(pts2) = c(x1, y1)
                

		
		if (plotMap==TRUE) {
			if(is.null(spMap)) stop(paste("if plotMap=TRUE, spMap must also be provided"))
			plot(spMap, xlim=c(bbox(spMap)[1],bbox(spMap)[3]), ylim=c(bbox(spMap)[2],bbox(spMap)[4]), ...)
		    plot(pts2, add=TRUE)
		}
		
g <- gstat(formula=formula, locations= loc, data=data, model = fitmodel, ...)
interp <- predict(g, newdata = pts2, BLUE = BLUE)
returnVal <- sum(sqrt(interp[["var1.var"]]))/n
returnVal
	}
	results <- rbga(as.matrix(c(rep(xmin,n), rep(ymin, n))), as.matrix(c(rep(xmax,n), rep(ymax, n))), popSize=popSize, evalFunc=evaluate, verbose=TRUE, iters=generations, ...)
return(results)
}
