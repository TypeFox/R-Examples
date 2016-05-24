plotBclim <-
function(x,dim=1,title=NULL,presentleft=TRUE,blob=TRUE,MDPcol="blue",denscol="red",MDPtransp=0.1,denstransp=0.5,leg=TRUE,legloc="topleft",...) {
  
  if(class(x)!="Bclim") stop("Needs a Bclim output object")  
  # Set up plot 
  par(mar=c(4,4,3,1))
  xrange <- range(c(0,x$time.grid))
  if(!presentleft) xrange <- rev(xrange)
  yrange <- range(c(0,x$MDP[,,dim]))
  mytitle <- title
  if(is.null(title)) mytitle <- paste(x$core.name,": ",x$clim.dims[dim],sep="")
  
  if(dim==1) {
    plot(1,1,type="n",xlim=xrange,ylim=yrange,xlab="Age (k cal years BP)",ylab=expression(paste("GDD5 (",degree,"C days)",sep="")),las=1,bty="n",main=mytitle,xaxt='n')
  }
  if(dim==2) {
    plot(1,1,type="n",xlim=xrange,ylim=yrange,xlab="Age (k cal years BP)",ylab=expression(paste("MTCO (",degree,"C)",sep="")),las=1,bty="n",main=mytitle,xaxt='n')    
  }
  if(dim==3) {
    plot(1,1,type="n",xlim=xrange,ylim=yrange,xlab="Age (k cal years BP)",ylab=x$clim.dims[3],las=1,bty="n",main=mytitle,xaxt='n')
  }

    axis(side=1,at=pretty(x$time.grid,n=10))
    rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col="lightgray",border="NA")
    grid(col="white")
    
    # Plot MDPs if required
    if(blob==TRUE) {
      # Find the one with the smallest amount of simulations
      num <- min(x$nchron,x$n.samp)
      # Get simulations of dates from chronology
      chron <- read.table(x$Chronsfile,nrows=num)  
      #chron <- x$chron.store
      
      # Sort out colours
      tmp <- col2rgb(MDPcol)
      mycol <- rgb(tmp[1,1]/255,tmp[2,1]/255,tmp[3,1]/255)
      mycol2 <- paste(mycol,as.character(as.hexmode(round(MDPtransp*255,0))),sep="")
              
      for(i in 1:x$n) {
        # Get simulations of MDP for dimension dim and layer i
        MDPtemp <- x$MDP[1:num,i,dim]
        chrontemp <- chron[1:num,i]
          
        if(var(chrontemp)>0) {
          tmp <- MASS::kde2d(chrontemp,MDPtemp)
          
          # Standardise and find 95% limit
          z <- tmp$z/sum(tmp$z)
          prop <- 1
          limitsdiff <- (max(z)-min(z))/100
          limits <- min(z)
          
          while(prop>0.95) {
            limits <- limits+limitsdiff
            prop <- sum(z[z>limits])
          }
          
          tmp2 <- contourLines(tmp$x,tmp$y,z,levels=limits)
          for(j in 1:length(tmp2)) polygon(tmp2[[j]]$x,tmp2[[j]]$y,col=mycol2,border=mycol2)  
        } else {
          lines(c(chrontemp[1],chrontemp[2]),quantile(MDPtemp,probs=c(0.025,0.975)),col=mycol2,lwd=3)
        }
      }
    }
    
    MAP <- rep(NA,length(x$time.grid))
    HDR95 <- matrix(NA,nrow=length(x$time.grid),ncol=50)
    HDR75 <- matrix(NA,nrow=length(x$time.grid),ncol=50)
    HDR50 <- matrix(NA,nrow=length(x$time.grid),ncol=50)
    for(i in 1:length(x$time.grid)) {
        if(sd(x$clim.interp[,i,dim])>0) {
          TempHDR <- hdrcde::hdr(x$clim.interp[,i,dim],h=bw.nrd0(x$clim.interp[,i,dim]),prob=c(1,50,75,95))$hdr
          
          Temp95 <- TempHDR[1,!is.na(TempHDR[1,])]
          HDR95[i,1:length(Temp95)] <- Temp95
          Temp75 <- TempHDR[2,!is.na(TempHDR[2,])]
          HDR75[i,1:length(Temp75)] <- Temp75
          Temp50 <- TempHDR[3,!is.na(TempHDR[3,])]
          HDR50[i,1:length(Temp50)] <- Temp50
          MAP[i] = mean(TempHDR[4,1:2])
        } else {
          HDR95[i,1:2] <- x$clim.interp[1,i,dim]
          HDR75[i,1:2] <- x$clim.interp[1,i,dim]
          HDR50[i,1:2] <- x$clim.interp[1,i,dim]
          MAP[i] = x$clim.interp[1,i,dim]
        }
    }

    
    if(blob==TRUE) {    
      # Sort out colours
      tmp <- col2rgb(denscol)
      mycol <- rgb(tmp[1,1]/255,tmp[2,1]/255,tmp[3,1]/255)
      mycol2 <- paste(mycol,as.character(as.hexmode(round(denstransp*255,0))),sep="")
    
      # Now draw a pretty polygon
      polygon(c(x$time.grid,rev(x$time.grid)),c(apply(HDR95,1,"min",na.rm=TRUE),rev(apply(HDR95,1,"max",na.rm=TRUE))),col=mycol2,border=mycol2)
      polygon(c(x$time.grid,rev(x$time.grid)),c(apply(HDR75,1,"min",na.rm=TRUE),rev(apply(HDR75,1,"max",na.rm=TRUE))),col=mycol2,border=mycol2)
      polygon(c(x$time.grid,rev(x$time.grid)),c(apply(HDR50,1,"min",na.rm=TRUE),rev(apply(HDR50,1,"max",na.rm=TRUE))),col=mycol2,border=mycol2)
        
    } else {
    
        for(i in 1:length(x$time.grid)) {
          # Draw some pretty lines
          for(j in 1:(ncol(HDR95)/2)) {
            lines(c(x$time.grid[i],x$time.grid[i]),c(HDR95[i,2*j-1],HDR95[i,2*j]),col="blue",lwd=1)
          }
          for(j in 1:(ncol(HDR75)/2)) {
            lines(c(x$time.grid[i],x$time.grid[i]),c(HDR75[i,2*j-1],HDR75[i,2*j]),col="blue",lwd=2)
          }
          for(j in 1:(ncol(HDR50)/2)) {
            lines(c(x$time.grid[i],x$time.grid[i]),c(HDR50[i,2*j-1],HDR50[i,2*j]),col="blue",lwd=3)
          }
        }
    }
    
    # Finally draw a legend
    if(leg==TRUE) {
      legend(legloc,legend=c("95/75/50% Joint Posteriors","95% Marginal Data Posteriors"),fill=c(denscol,MDPcol),bty="n")
    }
    
# End of function   
}
