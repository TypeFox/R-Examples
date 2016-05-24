plot.sizepower <- function(x,xnull,legend.pos=NULL,...) {

	Fx <- x
	
	Fxnull <- xnull
	
	if(!inherits(Fx, "Fx") | !inherits(Fxnull, "Fx")) stop("Method is only for Fx objects!")

	## Retrieve informations from Fx computed by calcFx()
	Fxi <- Fx$Fx.mat
	xi <- Fx$x
	law <- Fx$law
	statnames <- Fx$statnames
	N <- Fx$N
	
	## Retrieve informations from Fxnull computed by calcFx()
	Fxnull <- Fxnull$Fx

	## now we plot these (Fxnull,Fxi) in the same graph. Fxnull = null distribution, Fxi = alternative distribution

	# we trace the 45 degree line
	plot(c(0,1),c(0,1),type = "n",main=paste("Size-power curves","\n N = ",N," & m = ",length(xi),sep=""),
             sub=paste("Data: ",law,sep=""),xlab="size",ylab="power",asp = 1)
	lines(x=c(0,1),y=c(0,1),lty="dotted",col=1)
	text(0.9,0.85,"45 degree line",cex=1,srt=45)
	
	# then add theses Size-power curves
	for (i in 1:length(statnames)) {
          points(c(0,Fxnull[,i],1),c(0,Fxi[,i],1),type="l",col=i+1)
	}
	
	# legend
        if (is.character(legend.pos)) legend(x=legend.pos,legend=statnames,col=2:(length(statnames)+1),lty=1)
        if (is.null(legend.pos)) {
          legend.x <- 0.7
          legend.y <- 0.3
          legend(x=legend.x,y=legend.y,legend=statnames,col=2:(length(statnames)+1),lty=1)
        } else {
          legend(x=legend.pos[1],y=legend.pos[2],legend=statnames,col=2:(length(statnames)+1),lty=1)
        }

}
