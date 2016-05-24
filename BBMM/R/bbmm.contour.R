bbmm.contour <-
function(x, levels, locations=NULL, plot=TRUE){
    
	#
	# x = object of class bbmm
    # levels = vector of desired contour levels. E.g., c(85, 90, 95, 99)
    # locations = (optional, used if plot=TRUE) x and y coordinates of locations used to estimate Brownian bridge movement model
    # plot = logical
    #
	
  v <- na.omit(x$probability)
  probability <- NULL
  for (i in levels) {
		contour.z <- function(z) {
            abs(i/100 - sum(v[v >= z])/sum(v))
        }
    probability <- c(probability, optimize(contour.z, 
		c(0, max(v)), tol = .Machine$double.eps)$minimum)
  }

   ans <- vector("list", 2)
   names(ans) <- c("Contour", "Z")
   ans[[1]] <- paste(levels, "%", sep='')
   ans[[2]] <- probability
   if(plot){
        unique.xy <- expand.grid(x=unique(x[[2]]), y=unique(x[[3]]))
        temp <- data.frame(x=x[[2]], y=x[[3]], z=x[[4]])
        temp <- merge(unique.xy, temp, all.x=TRUE)
		temp <- temp[order(temp$x, temp$y),]
        temp$z[is.na(temp$z) == TRUE] <- 0 
        z <- matrix(temp$z, nrow=length(unique(temp$x)), 
					ncol=length(unique(temp$y)), byrow=TRUE)
        contour(x=unique(temp$x), y=unique(temp$y), z, levels=ans[[2]], 
                drawlabels=TRUE, labels=ans[[1]], xlab="X", ylab="Y")
        points(locations[,1], locations[,2], pch=19, cex=0.7, col="red")
        lines(locations[,1], locations[,2], col="red")
   }
   return(ans)
}
