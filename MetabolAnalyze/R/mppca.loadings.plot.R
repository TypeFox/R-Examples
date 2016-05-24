mppca.loadings.plot <-
function(output, Y, barplot=FALSE, labelsize=0.3)
{
  q<-output$q
  g<-output$g
  if(q==1)
  {
  	for(i in 1:g)
  	{
     barplot2(output$loadings[,,i], ylim=c(min(output$loadings[,,i])-1, max(output$loadings[,,i])+1), ylab="Loadings", xlab="Spectral regions",
     las=2, col="red", font.main = 2, plot.grid=TRUE, names.arg = colnames(Y), cex.names=labelsize, main=paste("Group ", i, sep=""))
     if(i < g)
     {
        ask(msg = "Press <RETURN> to view the next loadings plot: ")
     }
    }#i
  }
  
  if(q==2)
  {
     if(barplot==FALSE)
     {
     	for(i in 1:g)
     	{
         plot(output$loadings[,,i][,1], output$loadings[,,i][,2],xlab="PC 1", ylab="PC 2", type="n", main=paste("Group ", i, sep=""))
         text(output$loadings[,,i][,1], output$loadings[,,i][,2], colnames(Y), cex=1)
         abline(h=0, col="red", lty=1); abline(v=0, col="red")
         if(i < g)
     	 {
        	ask(msg = "Press <RETURN> to view the next loadings plot: ")
     	 }
        }#i
     }#if
     else{
     	for(k in 1:g)
     	{
     	 for(i in 1:q)
     	 {
          barplot2(output$loadings[,,k][,i], ylim=c(min(output$loadings[,,k][,i])-1,max(output$loadings[,,k][,i])+1), las=2, width=0.5, space=0.5, col=i, ylab="Loadings", xlab="Spectral regions", names.arg = colnames(Y), cex.names=labelsize, plot.grid=TRUE, font.main = 2, main=paste("Loadings on PC", i,  " for group ", k, sep=""))
		  ask(msg = "Press <RETURN> to view the next loadings plot: ")          
         } #i
        } #k
     }#else
  }#if
  
  if(q>2)
  {
  	if(barplot==FALSE)
  	{
  		for(k in 1:g)
  		{
		 pairs(output$loadings[,,k], lower.panel=function(x,y,names){text(x,y,colnames(Y), cex=0.8);abline(h=0, col="red", lty=1); abline(v=0, col="red")}, upper.panel=function(x,y,names){text(x,y,colnames(Y), cex=0.8);abline(h=0, col="red", lty=1); abline(v=0, col="red")}, labels=paste("PC", 1:q, sep=""), main=paste("Loadings for group ", k, sep=""))
		 ask(msg = "Press <RETURN> to view the next loadings plot: ")
		 }
	}else{
	   for(k in 1:g)
	   {
		for(i in 1:q)
		{
	    	barplot2(output$loadings[,,k][,i], ylim=c(min(output$loadings[,,k][,i])-1,max(output$loadings[,,k][,i])+1), las=2, width=0.5, space=0.5, col=i, ylab="Loadings", xlab="Spectral regions", names.arg = colnames(Y), cex.names=labelsize, plot.grid=TRUE, font.main = 2, main=paste("Loadings on PC", i, " group", k, sep=""))
            ask(msg = "Press <RETURN> to view the next loadings plot: ")            
        } #i
       }#k
	} # else
  }
} # End plot.mppca.loadings

