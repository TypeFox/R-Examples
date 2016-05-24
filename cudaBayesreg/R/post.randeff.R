post.randeff <-
function(out, classnames=NULL, climits=TRUE)
{
		nvar <- dim(out$betadraw)[2]
		nreff <- ncol(out$Deltadraw)
		# Deltadraw: Posterior density of random effects
		ng <- nreff/nvar # total number of groups with intercept
		cc <- seq(1:nvar)
		i1=ng*c(0:(nvar-1))+1 # draws of the mean for all groups
		matplot(out$Deltadraw[,i1], type="l", col=cc, main=expression(
			paste("Draws of the mean of the random effects distribution,",Delta)))
		mtext(paste("Draws for each of the",nvar,"regression variables"),3)
		legend("topleft", paste("var",cc), col=cc, lty = c(1, 1, 1), merge = TRUE)
		# class plots  
		if(is.null(classnames)) return();
		nshow <- length(classnames)+1
    ylim <- c(NULL,NULL)
    if(climits)
      ylim <- range(out$Deltadraw[-i1])
		if( nreff <= nvar ) {
			cat("Unit random effects not available for Z=iota\n")
		}
		else {
			par(ask=T)
			for(k in 2:nshow) {
				gi=ng*c(0:(nvar-1))+k
				matplot(out$Deltadraw[,gi], type="l", col=cc, ylim=ylim,
				 main=paste("Random effects for class ",classnames[k-1],sep=""))
				mtext(paste("Draws for each var of the",nvar,"regression variables"),3)
				legend("topleft", paste("var",cc), col=cc, lty = c(1, 1, 1), merge = TRUE)
			}
			par(ask=F)
		}
}

