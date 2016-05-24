plot.lad <-
function(x,...)
{
	ans <- summary.lad(x);

	if (!ans$yfactor)
	{
		dd <- NCOL(ans$R);

		if (dd==1)
		{ 
			plot(ans$y, ans$R, pch=20, ylab="First Component", xlab="y", ...); 
			lines(lowess(ans$y, ans$R), lwd=2, col="red");
			title(main="Response against the reduction")
		}
		if (dd==2) 
		{
			par(mfrow=c(1,2), oma=c(0,0,2,0)); 
			plot(ans$y, ans$R[,1], pch=20, ylab="First Component", xlab="y",...); 
			lines(lowess(ans$y, ans$R[,1]), lwd=2, col="red");
			plot(ans$y, ans$R[,2], pch=20, ylab="Second Component", xlab="y",...);
			lines(lowess(ans$y, ans$R[,2]), lwd=2, col="red");
			title(main="Components of the reduction vs Response", outer=TRUE)		
		}
		if (dd==3) 
		{
			par(mfrow=c(2,2), oma=c(0,0,2,0)); 
			plot(ans$y, ans$R[,1], pch=20, ylab="First Component", xlab="y",...); 
			lines(lowess(ans$y, ans$R[,1]), lwd=2, col="red");

			plot(ans$y, ans$R[,2], pch=20, ylab="Second Component", xlab="y",...); 
			lines(lowess(ans$y, ans$R[,2]), lwd=2, col="red");

			plot(ans$y, ans$R[,3], pch=20, ylab="Third Component", xlab="y",...);
			lines(lowess(ans$y, ans$R[,3]), lwd=2, col="red");
			title(main="Components of the reduction vs Response", outer=TRUE)
		}
		if (dd>3) 
		{
			par(mfrow=c(2,2), oma=c(0,0,2,0)); 

			plot(ans$y, ans$R[,1], pch=20, ylab="First Component", xlab="y",...); 
			lines(lowess(ans$y, ans$R[,1]), lwd=2, col="red");

			plot(ans$y, ans$R[,2], pch=20, ylab="Second Component", xlab="y",...); 
			lines(lowess(ans$y, ans$R[,2]), lwd=2, col="red");

			plot(ans$y, ans$R[,3], pch=20, ylab="Third Component", xlab="y",...); 
			lines(lowess(ans$y, ans$R[,3]), lwd=2, col="red");

			plot(ans$y, ans$R[,4], pch=20, ylab="Fourth Component", xlab="y",...);
			lines(lowess(ans$y, ans$R[,4]), lwd=2, col="red");

			title(main="Components of the reduction vs Response", outer=TRUE)
		}
	}
	if (ans$yfactor)
	{
		mycolors <- c("black", "blue", "red", "green", "yellow", "gray", "cyan", "magenta");
		pchs <- c(20, 21, 22, 17, 1, 2, 3, 4);
		nl <- length(unique(ans$y));
		mycol <- mypch <- as.integer(factor(ans$y, levels=unique(ans$y)));
		for (i in 1:nl)  { mycol[mycol==i]<- mycolors[i]; mypch[mypch==i] <- pchs[i]}

		if (NCOL(ans$R)==1)
		{
			par(mfrow=c(1,2), oma=c(0,0,2,0)); 
			plot(x=c(-1,1), xlab="", ylab="", ylim=c(0,1), xlim=c(0,1), col="white", cex=1.5, xaxt="n", yaxt="n", axes=FALSE);
			legend(0.2, 1, unique(ans$y), border = "blank", title = "Legend", pch=pchs[1:nl], col=mycolors[1:nl])

			plot(ans$R[,1], xlab="index", ylab="LAD - Dir1", pch=mypch, col=mycol, cex=1.5);
			title(main="Scatterplot of the sufficient reduction", outer=TRUE)
		}
		if (NCOL(ans$R)==2)
		{
			par(mfrow=c(1,2), oma=c(0,0,2,0)); 
			plot(x=c(-1,1), xlab="", ylab="", ylim=c(0,1), xlim=c(0,1), col="white", cex=1.5, xaxt="n", yaxt="n", axes=FALSE);
			legend(0.2, 1, unique(ans$y), border = "blank", title = "Legend", pch=pchs[1:nl], col=mycolors[1:nl])

			plot(ans$R[,1], ans$R[,2], xlab="LAD - Dir1", ylab="LAD - Dir2", pch=mypch, col=mycol)
			title(main="Scatterplot of the components of the sufficient reduction", outer=TRUE)
		}
		if (NCOL(ans$R) == 3)
		{
			par(mfrow=c(2,2), oma=c(0,0,2,0)); 
			plot(x=c(-1,1), xlab="", ylab="", ylim=c(0,1), xlim=c(0,1), col="white", cex=1, xaxt="n", yaxt="n", axes=FALSE);
			legend(0.2, 1, unique(ans$y), border = "blank", cex = 0.9, title = "Legend", pch=pchs[1:nl], col=mycolors[1:nl])

			plot(ans$R[,1], ans$R[,2], xlab="LAD - Dir1", ylab="LAD - Dir2", pch=mypch, col=mycol);
			plot(ans$R[,1], ans$R[,3], xlab="LAD - Dir1", ylab="LAD - Dir3", pch=mypch, col=mycol);
			plot(ans$R[,2], ans$R[,3], xlab="LAD - Dir2", ylab="LAD - Dir3", pch=mypch, col=mycol);
			title(main="Scatterplots of the components of the sufficient reduction", outer=TRUE)
		}
		if (NCOL(ans$R) > 3)
		{
			par(mfrow=c(3,3), oma=c(0,0,2,0));  
			plot(x=c(-1,1), xlab="", ylab="", ylim=c(0,1), xlim=c(0,1), col="white", cex=0.8, xaxt="n", yaxt="n", axes=FALSE);
			legend(0.2, 1, unique(ans$y), border = "blank", cex = 0.7, title = "Legend", pch=pchs[1:nl], col=mycolors[1:nl])

			plot(ans$R[,1], ans$R[,2], xlab="LAD - Dir1", ylab="LAD - Dir2", pch=mypch, col=mycol);
			plot(ans$R[,1], ans$R[,3], xlab="LAD - Dir1", ylab="LAD - Dir3", pch=mypch, col=mycol);
			plot(ans$R[,1], ans$R[,4], xlab="LAD - Dir2", ylab="LAD - Dir3", pch=mypch, col=mycol);
			plot(ans$R[,2], ans$R[,3], xlab="LAD - Dir2", ylab="LAD - Dir3", pch=mypch, col=mycol);
			plot(ans$R[,2], ans$R[,4], xlab="LAD - Dir2", ylab="LAD - Dir3", pch=mypch, col=mycol);
			title(main="Scatterplots of the components of the sufficient reduction", outer=TRUE)
		}
	}
}
