plot.pvalue <- function(x,legend.pos=NULL,...) {

	Fx <- x
	
	if(!inherits(Fx, "Fx")) stop("Method is only for Fx objects!")

	## Retrieve informations from Fx computed by calcFx()
	Fxi <- Fx$Fx.mat
	xi <- Fx$x
	law <- Fx$law
	statnames <- Fx$statnames
	N <- Fx$N

	## Now we plot these (xi,Fxi) in the same graph

	# we trace the 45 degree line
	plot(c(0,1),c(0,1),type = "n",main=paste("P value plots","\n N = ",N," & m = ",length(xi),sep=""),
             sub=paste("Data: ",law,sep=""),xlab="Nominal Size",ylab="Actual Size",asp = 1,...)
	lines(x=c(0,1),y=c(0,1),lty="dotted",col=1)
	text(0.9,0.85,"45 degree line",cex=1,srt=45)
	
	# then add theses P value plots
	for (i in 1:length(statnames)) {
          points(xi,Fxi[,i],type="l",col=i+1,...)
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
