.plot.EPSGO.parms<- function (Xtrain, Ytrain,bounds, Ytrain.exclude=10^16, plot.name=NULL ){

# for interval search plot visited points and the Q values (=Ytrain)
# exclude: for D=1 make  an additional plot: skip values for empty model, for example: Ytrain.exclude=10^16
	
	
	if (is.null(plot.name)) plot.name<- "Latin hypercube sampling"
	
	# dimension of interval
	D<-ncol(Xtrain)
	
	if (D> 2) print("Warning! Plot only one or two dimension parameter spaces.") 
		if ( D<=2 ){
		# 1D plot
		if (D==1) {
			par(mfrow=c(1,3))
			plot(Xtrain, xlab="Index", ylab=rownames(bounds)[1], col="orange", pch=20,
						 main=paste( plot.name,", n=",nrow(Xtrain)*D) )
			abline(h=seq(bounds[1,1],bounds[1,2],length=(nrow(Xtrain)+1)), lty=2, col=3)
		} else {
			# 2D plot
			plot(Xtrain, xlab=rownames(bounds)[1], ylab=rownames(bounds)[2], col="orange", pch=20, 
						 main=paste( plot.name,", n=",nrow(Xtrain) ) )
			abline(v=seq(bounds[1,1],bounds[1,2],length=(nrow(Xtrain)+1)), lty=2, col=3)
			abline(h=seq(bounds[2,1],bounds[2,2],length=(nrow(Xtrain)+1)), lty=2, col=3)
		}
	}
	# add start.q.values to the plot
	if ( D ==1)  text( c(1: nrow(Xtrain)),Xtrain, labels=round(Ytrain,6), pos=1, cex=0.5 )
	if ( D ==2)  text( Xtrain, labels=signif(Ytrain,6), pos=1, cex=0.5 )
	
	# for D=1  Ytrain ~ Xtrain plot
	
	if ( D ==1) {
		# sort X and Y to have a cont. line
		tmp<- cbind(Xtrain,Ytrain)	
		tmp.sort<- sortmat(tmp,1, decreasing=FALSE)
		plot(tmp.sort, type="b")
		
		if (!is.null(Ytrain.exclude)){
			# if exclude for Ytrain exists, skip those  points 
			tmp.sort2<- data.frame(tmp.sort)
			tmp.sort2$Ytrain[tmp.sort2$Ytrain == Ytrain.exclude]<- NA
			
			plot(tmp.sort2, type="b", main=paste("excluded Ytrain=", Ytrain.exclude ) )
		}
		
	}

}

	
	