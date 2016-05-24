

plot.ir_curve <- function(x,ylim=c(),xlim=c(),lwd=2, type="l",
				xlab="Maturity (years)",ylab="Zero-coupon yields (in percent)", 
				col="steelblue",lty=1, ...) 
				{
	plot(x[,1] ,x[,2]*100, type=type, ylim=ylim, xlim=xlim, xlab=xlab,
     ylab=ylab,lwd=lwd,lty=lty,col=col, ... )
      
}


plot.spot_curves <- function(x,multiple= FALSE,
				ylim= c(range(mapply(function(i) 
				range(x[[i]][,2]),seq(x))))*100,xlim=c(),
				type="l", lty=1, lwd=2, expoints=NULL, 
				ylab= "Zero-coupon yields (percent)",
				xlab= "Maturity (years)",main="Zero-coupon yield curves",
				...) {

	if(multiple) 
	{ plot(x[[which.max(mapply(function(i) max(x[[i]][,1]),
		seq(x)))]][,1], x[[which.max(mapply(function(i) 
		max(x[[i]][,1]), seq(x)))]][,2]*100, 
        type=type,col=which.max(mapply(function(i) max(x[[i]][,1]),seq(x))),
        lty=lty,lwd=lwd,xlab=xlab,ylab=ylab,ylim=ylim, ... )
   
	  for(k in c((seq(x))[-which.max(mapply(function(i) max(x[[i]][,1]), seq(x)))]))
	  { lines(x[[k]][(if(is.numeric(expoints)) seq(expoints[k]) 
	  	else seq(nrow(x[[k]]))),1],x[[k]][(if(is.numeric(expoints)) seq(expoints[k]) 
	  	else seq(nrow(x[[k]]))),2]*100,col=k,lwd=lwd,lty=lty, ... )
      
        if(is.numeric(expoints))
        {lines(x[[k]][((expoints[k]+1):nrow(x[[k]])) ,1],
        	x[[k]][((expoints[k]+1):nrow(x[[k]])),2]*100,col=k,lwd=lwd,lty=5, ... )
         }
    title(main)
    legend("bottom",legend=names(x),col=seq(x),lty=lty,lwd=lwd)
   }
  }
	else
	{   old.par <- par(no.readonly = TRUE)
		par(ask=TRUE)
		for(k in seq(x)) 
		{ 
		plot.ir_curve(x[[k]],...)
		title(names(x)[k])
      	legend("bottom",legend=main,col=c("steelblue"), lty = 1 , pch=c(-1))
        }	
        on.exit(par(old.par))	
	}
	

}


plot.fwr_curves <- function(x,multiple= FALSE,
					ylim= c(range(mapply(function(i) range(x[[i]][,2]),
					seq(x))))*100,xlim=c(),type="l", lty=1, 
					lwd=2, expoints=NULL, ylab= "Forward rate (percent)",
					xlab= "Maturity (years)",main="Forward rate curves",...) 
					
{ plot.spot_curves(x,ylab=ylab, xlab=xlab, main=main,
	multiple=multiple,expoints=expoints,lty=lty,lwd=lwd,type=type, ... )

}


plot.s_curves <- function(x,xlim=c(range(mapply(function(i) 
					range(x[[i]][,1]),seq(x)))),
					ylim=c(range(mapply(function(i) range(x[[i]][,2]),
					seq(x))))*10000,expoints=NULL, xlab="Maturity (years)", 
					ylab="Spread (basis points)", lwd=2,lty=1, main="Spread curves", ...)
					
{  if(!is.character(x))
   {
	
   plot(0,0, type="n",xlab=xlab,ylab=ylab,xlim=xlim, ylim=ylim,...)
 
   for(k in c(2:length(x)))
    { lines(x[[k]][(if(is.numeric(expoints))
    	seq(expoints[k]) else seq(nrow(x[[k]]))),1],
    	x[[k]][(if(is.numeric(expoints)) seq(expoints[k]) 
    	else seq(nrow(x[[k]]))),2]*10000,col=k,lwd=lwd,lty=lty, ... )
      
      if(is.numeric(expoints))
      {lines(x[[k]][((expoints[k]+1):nrow(x[[k]])) ,1],
     	x[[k]][((expoints[k]+1):nrow(x[[k]])),2]*10000,
     	col=k,lwd=lwd,lty=5, ... )
   
       }
    } 
   title(main)
   legend("topleft",legend=names(x)[-1],col=seq(x)[-1],lty=1,lwd=lwd)
   } else warning("No spread curves available")
} 


plot.error <- function(x,type="b",main="", mar= c(7,6,6,2) + 0.1, oma=c(4,2,2,2) +0.1,
                       ylab="Error", ...) {
	old.par <- par(no.readonly = TRUE)
    par(mar=mar, oma=oma, ... )
    
		plot(x[,1],x[,2],axes=FALSE,pch=19,lwd=c(1,2),xlab="", ylab=ylab,type=type, ...)
		axis(1,x[,1],rownames(x),las=3,...)
		axis(2,...)
		axis(3,x[,1],round(x[,1],2),...)
		mtext("Maturity (years)",3,line=2.5)
		lines(x[,1],rep(0,nrow(x)),lty=2,lwd=1,... )
		title(xlab="ISIN", outer=TRUE,main=main,...) 
	 
	 on.exit(par(old.par))
}   


plot.df_curves <- function(x,multiple= FALSE,
					ylim= c(range(mapply(function(i) range(x[[i]][,2]),
					seq(x))))*100,xlim=c(),type="l", lty=1,
					lwd=2, expoints=NULL, ylab="Discount factor (percent)",
					xlab= "Maturity (years)",main="Discount factor curves",...) 
{ plot.spot_curves(x,ylab=ylab, xlab=xlab, main=main,
	multiple=multiple,expoints=expoints,lty=lty,lwd=lwd,type=type, ... )

		}
