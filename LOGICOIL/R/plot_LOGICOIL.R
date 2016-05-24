plot_LOGICOIL <-
function(prob.oligo, id)
{

	png(sprintf(filename="%s_LOGICOIL.png", id), height=800, width=1200, bg="white", res=300)
	par(mar=c(3.5,3.5,0.75,0.75),family="Arial", font=1, cex.lab=0.75, cex.axis=0.75) 	  	
	mymax <- max(as.numeric(prob.oligo)-1) + 0.1
	pal=c("#304E67","darkgreen","tomato","blue")	
	
	plot(as.numeric(prob.oligo[,1])-1,
      type="h",
      col=pal[1],
	    ylim=c(-mymax, mymax),
      lwd=1.5,bty="l",
	    axes=FALSE,
      ann=FALSE,
      cex=1.0,
      tck=1
	)
	lines(as.numeric(prob.oligo[,2])-1, type="h",  col=pal[2], lwd=1.5)
	lines(as.numeric(prob.oligo[,3])-1, type="h",  col=pal[4], lwd=1.5)
	lines(as.numeric(prob.oligo[,4])-1, type="h", col=pal[3], lwd=1.5)
	xlimit <- c()
	xlimit <- ceiling(nrow(prob.oligo)/5)*5
	axis(1, lwd.ticks=0, at=seq(0,xlimit,5), labels=seq(0,xlimit,5))
	axis(2, lwd=0, las=1, at=seq(-1,1,0.25), yaxp=c(0,100,4))
	
	for(y in seq(-1,1,0.25))
	{
	  lines(rep(y, nrow(prob.oligo)), type="s", col="lightgray", lwd=1)
	}
	
	mtext("Sequence id",side=1,line=2.5,cex=0.5)
	mtext("LOGICOIL scores",side=2,line=2.5,cex=0.5)
	legend("topright",ncol=2, c("Antiparallel Dimer","Parallel Dimer","Trimer","Higher-order"), cex=0.5, 
	    col=pal, 
	    lty=1:1, # solid 
	    lwd=2, # line width
	    bty="n"
	)
	dev.off()
}
