             
overlayCI <- function(cis, xpos=NULL,ci.cols=NULL,  ci.ex=2, ci.ocol = "grey40", p.col="grey40",pch=1, sig.col = "red",sig.lwd = 1,yusr=NULL,ci.label="Differences",ci.cex=0.5,arrow.length=0.1,...)
{    
	 # draw the confidence intervals
      npairs <- nrow(cis)
      if (is.null(xpos))
          xp <- (1:npairs)+.5
      else xp <- xpos
      if (is.null(yusr)) {
        yrange <- pretty(range(cis,na.rm=TRUE))
        yusr <- c(min(yrange), max(yrange))
        }
      xusr <- par("usr")[1:2]
      usr <- c(xusr,yusr)
        
      par("usr"=usr)
      
      m <- cis[,1]
      cis <- cis[,-1]
	  nlevels <- ncol(cis)/2
      if (is.null(ci.cols)) ci.cols <- gray.colors(nlevels, start=0.4, end=0.8)
      n <- ncol(cis) 
      cis <- cis[,n:1]
      for (lev in 1:nlevels) {
      	  j <- 2*(lev-1)
      	  
	      for (i in 1:length(xp)) 
	      if (!is.na(cis[i,1])) {
	              lines(c(xp[i],xp[i]),
	                  cis[i,(1:2) + j],
	                  col=ci.cols[lev],
	                  lwd= ci.ex * lev,lend=2)
		      	   }
	      }
        for (i in 1:length(xp)) 
         if (!is.na(cis[i,1]))
       	  lines(c(xp[i]-.2, xp[i]+.2), c(0,0),col=ci.ocol)
       	if (!is.null(ci.ocol)) arrows(usr[2],0, usr[2]-.5,0, col=ci.ocol,length=arrow.length)

       	points(xp, m, col=p.col,pch=pch,cex=ci.cex)       	      
       	 limits <- cis[,2:1]  
       
          for (i in 1:npairs)	{
          
          if (!is.na(limits[i,1]) && limits[i,1]>=0) {
		       arrows(xp[i],-limits[i,1], xp[i],0,
	           col=sig.col, length=arrow.length,lwd=sig.lwd,
	          lend=2,lty=1)
		  }
		  if (!is.na(limits[i,1]) && limits[i,2]<=0) {
	       arrows(xp[i], -limits[i,2],xp[i], 0,
	          col=sig.col, length=arrow.length,lwd=sig.lwd,
	          lend=2,lty=1)
	      }}
	     # arrows(length(o)+1,0, length(o)+.5,0, col=ci.ocol,lwd=2)
          axis(4,col=ci.ocol,col.axis=ci.ocol)
          axis(4,labels=ci.label,at=0,line=par("mgp")[1]-1,tcl=0, cex.axis=par("cex.lab"),col=ci.ocol,col.axis=ci.ocol)
          invisible()
	      }
	      
	      

mc_plot <- function(data, fit,path =eulerian, col=rainbow(length(data),s=.4), levels=c(0.90,0.95,0.99), 
cifunction=function(a,level) TukeyHSD(a,conf.level= level)[[1]],
varwidth=TRUE,frame.plot=FALSE,boxwex=.3,cex=0.75,zoom=NULL, ci.yusr=NULL,ci.pos=FALSE,...)  {

	if (is.function(path)){
		tuk <- cifunction(fit, 0.95)
      # put p-values into distance matrix and order them
      d <- edge2dist(tuk[,4])
      o <- path(d)
      }
      else o <- path
      if (is.null(o)) o <- 1:length(data)
      if (!is.null(zoom)) {
    	o <- o[zoom]
    	}
       bp <-boxplot(data[o],col=col[o],varwidth=varwidth,
             frame.plot=frame.plot,boxwex=boxwex,cex=cex,...)
      
      nlevels <- length(levels)
      cis <- NULL
      for (lev in levels) {
	      tuk <- cifunction(fit,lev)
	      cis <- cbind(cis,tuk[,2:3])
	  }
	  if (!is.null(cis)){
      cis <-cbind(tuk[,1],cis)
      cis <- path_cis(cis,o,ci.pos=ci.pos)
      overlayCI(cis,yusr=ci.yusr,...) 
     }
      invisible()    		
	}	      
