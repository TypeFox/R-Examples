plot.pfc <-
function(x,...)
{
	if (is.numeric(x$y))
	{
		dd <- NCOL(x$R);

		if (dd==1)
		{ 
			plot(x$y, x$R, pch=20, ylab="First Component", xlab="y", ...); 
			lines(lowess(x$y, x$R), lwd=2, col="red");
	
			title(main="Response against the reduction")
		}

		if (dd==2) 
		{
			par(mfrow=c(1,2), oma=c(0,0,2,0)); 

			plot(x$y, x$R[,1], pch=20, ylab="First Component", xlab="y",...); 
			lines(lowess(x$y, x$R[,1]), lwd=2, col="red");

			plot(x$y, x$R[,2], pch=20, ylab="Second Component", xlab="y",...);
			lines(lowess(x$y, x$R[,2]), lwd=2, col="red");

			title(main="Components of the reduction vs Response", outer=TRUE)		
		}

		if (dd==3) 
		{
			par(mfrow=c(2,2), oma=c(0,0,2,0)); 

			plot(x$y, x$R[,1], pch=20, ylab="First Component", xlab="y",...); 
			lines(lowess(x$y, x$R[,1]), lwd=2, col="red");

			plot(x$y, x$R[,2], pch=20, ylab="Second Component", xlab="y",...); 
			lines(lowess(x$y, x$R[,2]), lwd=2, col="red");

			plot(x$y, x$R[,3], pch=20, ylab="Third Component", xlab="y",...);
			lines(lowess(x$y, x$R[,3]), lwd=2, col="red");

			title(main="Components of the reduction vs Response", outer=TRUE)
		}

		if (dd>3) 
		{
			par(mfrow=c(2,2), oma=c(0,0,2,0)); 

			plot(x$y, x$R[,1], pch=20, ylab="First Component", xlab="y",...); 
			lines(lowess(x$y, x$R[,1]), lwd=2, col="red");

			plot(x$y, x$R[,2], pch=20, ylab="Second Component", xlab="y",...); 
			lines(lowess(x$y, x$R[,2]), lwd=2, col="red");

			plot(x$y, x$R[,3], pch=20, ylab="Third Component", xlab="y",...); 
			lines(lowess(x$y, x$R[,3]), lwd=2, col="red");

			plot(x$y, x$R[,4], pch=20, ylab="Fourth Component", xlab="y",...);
			lines(lowess(x$y, x$R[,4]), lwd=2, col="red");

			title(main="Components of the reduction vs Response", outer=TRUE)
		}
	}

	if (is.factor(x$y))
	{
		mycolors <- c("black", "blue", "red", "green", "yellow", "gray", "cyan", "magenta")
		pchs <- c(20, 21, 22, 17, 1, 2, 3, 4)
		nl <- length(unique(x$y))
		mycol <- mypch <- as.integer(factor(x$y, levels=unique(x$y)))
		for (i in 1:nl)  { mycol[mycol==i]<- mycolors[i]; mypch[mypch==i] <- pchs[i]}

		if (NCOL(x$R)==1)
		{
			par(mfrow=c(1,2), oma=c(0,0,2,0)) 
			plot(x$R[,1], xlab="", ylab="", ylim=c(0,1), xlim=c(0,1), col="white", cex=1.5, xaxt="n", yaxt="n", axes=FALSE)
			legend(0.2, 0.8, unique(x$y), border = "blank", cex = 1, title = "Legend", pch=pchs[1:nl], col=mycolors[1:nl])

			plot(x$R[,1], xlab="index", ylab="PFC - Dir1", pch=mypch, col=mycol, cex=1.5)
			title(main="Scatterplot of the sufficient reduction", outer=TRUE)
		}


		if (NCOL(x$R)==2)
		{
			par(mfrow=c(1,2), oma=c(0,0,2,0)) 
			plot(x$R[,1], xlab="", ylab="", ylim=c(0,1), xlim=c(0,1), col="white", cex=1.5, xaxt="n", yaxt="n", axes=FALSE)
			legend(0.2, 0.8, unique(x$y), border = "blank", cex = 1, title = "Legend", pch=pchs[1:nl], col=mycolors[1:nl])

			plot(x$R[,1], x$R[,2], xlab="PFC - Dir1", ylab="PFC - Dir2", pch=mypch, col=mycol, cex=1)
			title(main="Scatterplot of the components of the sufficient reduction", outer=TRUE)
		}

		if (NCOL(x$R) == 3)
		{
			par(mfrow=c(2,2), oma=c(0,0,2,0)); 
			plot(x$R[,1], xlab="", ylab="", ylim=c(0,1), xlim=c(0,1), col="white", cex=1.5, xaxt="n", yaxt="n", axes=FALSE)
			legend(0.2, 0.8, unique(x$y), border = "blank", cex = 1, title = "Legend", pch=pchs[1:nl], col=mycolors[1:nl])

			plot(x$R[,1], x$R[,2], xlab="PFC - Dir1", ylab="PFC - Dir2", pch=mypch, col=mycol, cex=1)
			plot(x$R[,1], x$R[,3], xlab="PFC - Dir1", ylab="PFC - Dir3", pch=mypch, col=mycol, cex=1)
			plot(x$R[,2], x$R[,3], xlab="PFC - Dir2", ylab="PFC - Dir3", pch=mypch, col=mycol, cex=1)
			title(main="Scatterplots of the components of the sufficient reduction", outer=TRUE)
		}

		if (NCOL(x$R) > 3)
		{
			par(mfrow=c(3,3), oma=c(0,0,2,0));  
			plot(x$R[,1], xlab="", ylab="", ylim=c(0,1), xlim=c(0,1), col="white", cex=1.5, xaxt="n", yaxt="n", axes=FALSE)
			legend(0.2, 0.8, unique(x$y), border = "blank", cex = 1, title = "Legend", pch=pchs[1:nl], col=mycolors[1:nl])

			plot(x$R[,1], x$R[,2], xlab="PFC - Dir1", ylab="PFC - Dir2", pch=mypch, col=mycol, cex=1)
			plot(x$R[,1], x$R[,3], xlab="PFC - Dir1", ylab="PFC - Dir3", pch=mypch, col=mycol, cex=1)
			plot(x$R[,1], x$R[,4], xlab="PFC - Dir2", ylab="PFC - Dir3", pch=mypch, col=mycol, cex=1)
			plot(x$R[,2], x$R[,3], xlab="PFC - Dir2", ylab="PFC - Dir3", pch=mypch, col=mycol, cex=1)
			plot(x$R[,2], x$R[,4], xlab="PFC - Dir2", ylab="PFC - Dir3", pch=mypch, col=mycol, cex=1)
			title(main="Scatterplots of the components of the sufficient reduction", outer=TRUE)
		}

	}
}
