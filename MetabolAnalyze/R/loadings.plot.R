loadings.plot <-
function(output, barplot=FALSE, labelsize=0.3)
{
  q<-output$q
  if(q==1)
  {
     barplot2(output$loadings[,1], ylim=c(min(output$loadings)-1, max(output$loadings)+1), ylab="Loadings", xlab="Spectral regions",
     las=2, col="red", font.main = 2, plot.grid=TRUE, names.arg = rownames(output$loadings), cex.names=labelsize)
  }
  if(q==2)
  {
     if(barplot==FALSE)
     {
         plot(output$loadings[,1], output$loadings[,2],xlab="PC 1", ylab="PC 2", type="n", )
         text(output$loadings[,1], output$loadings[,2], rownames(output$loadings), cex=1)
         abline(h=0, col="red", lty=1); abline(v=0, col="red")
     }#if
     else{
     	for(i in 1:q)
     	{
          barplot2(output$loadings[,i], ylim=c(min(output$loadings[,i])-1,max(output$loadings[,i])+1), las=2, width=0.5, space=0.5, col=i, ylab="Loadings", xlab="Spectral regions", names.arg = rownames(output$loadings), cex.names=labelsize, plot.grid=TRUE, font.main = 2, main=paste("Loadings on PC", i, sep=""))
		   if((q > 1) & (i < q))
           {
              ask(msg = "Press <RETURN> to view the next loadings plot: ")
           }
        } #i
     }#else
  }#if
  
  if(q>2)
  {
  	if(barplot==FALSE)
  	{
		pairs(output$loadings, lower.panel=function(x,y,names){text(x,y,rownames(output$loadings), cex=0.8);abline(h=0, col="red", lty=1); abline(v=0, col="red")}, upper.panel=function(x,y,names){text(x,y,rownames(output$loadings), cex=0.8);abline(h=0, col="red", lty=1); abline(v=0, col="red")}, labels=paste("PC", 1:output$q, sep=""))
	}else{
		for(i in 1:q)
		{
			barplot2(output$loadings[,i], ylim=c(min(output$loadings[,i])-1,max(output$loadings[,i])+1), las=2, width=0.5, space=0.5, col=i, ylab="Loadings", xlab="Spectral regions", names.arg = rownames(output$loadings), cex.names=labelsize, plot.grid=TRUE, font.main = 2, main=paste("Loadings on PC", i, sep=""))
		   if((q > 1) & (i < q))
             {
              ask(msg = "Press <RETURN> to view the next loadings plot: ")
             }
        } #i
	} # else
  }
} # end plot.loadings

