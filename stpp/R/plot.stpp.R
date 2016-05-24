
plot.stpp <- function(x, s.region=NULL, t.region=NULL, mark=FALSE, mark.cexmin=0.4, mark.cexmax=1.2, mark.col=1, ...)
{
if (inherits(x,"stpp")==TRUE) 
	{ 
	if (mark==FALSE)
	  {
	  par(mfrow=c(1,2),pty="s")
	  if (is.null(s.region))	
	  plot(x[,1:2],main="xy-locations",...)
	  else
		{
		  polymap(s.region,xlab="x",ylab="y")
		  points(x[,1:2],...)
		  title("xy-locations")	 
		}
        xx=sort(x[,3],index.return=TRUE) 
        x=x[xx$ix,]
	  plot(x[,3],cumsum(x[,3]),type="l",xlab="t",ylab="",main="cumulative number",las=1,xlim=t.region)
 	  }
	if (mark==TRUE)
	 {
	  l=dim(x)[1]
	  CEX=seq(mark.cexmin,mark.cexmax,length=l)
        if (mark.col==0)
	     {
          par(mfrow=c(1,1),pty="s")
		  if (is.null(s.region))	
	  	   plot(x[,1:2],cex=CEX,...)
 		  else
			{
			  polymap(s.region,xlab="x",ylab="y")
			  points(x[,1:2],cex=CEX,...)	 
			}
          }
        else 
         {  
	       if (mark.col=="black" | mark.col==1)	
	           COL=grey((l:1)/l)
           if (mark.col=="red" | mark.col==2)	
	           COL=rgb(l:0, 0, 0, maxColorValue = l)
      	   if (mark.col=="green" | mark.col==3)	
        	   COL=rgb(0, l:0, 0, maxColorValue = l)
	       if (mark.col=="blue" | mark.col==4)	
	           COL=rgb(0, 0, l:0, maxColorValue = l)
	       par(mfrow=c(1,1),pty="s")
    	  if (is.null(s.region))	
	     	plot(x[,1:2],col=COL,cex=CEX,...)
 	      else
		  {
		   polymap(s.region,xlab="x",ylab="y")
		      points(x[,1:2],col=COL,cex=CEX,...)	 
		  }
         }
	 }	
      }
}
getS3method("plot", "stpp", optional = FALSE)

