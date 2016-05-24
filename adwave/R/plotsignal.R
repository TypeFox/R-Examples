plotsignal <-
function(x, ind= NULL, popA=NULL, popB=NULL, xlab=NULL, ylab=NULL, ylim=NULL, main=NULL){
  
  #check input arguments
  ctmp <- class(x)
  if (is.null(ctmp))  stop("object has no class")
  if (ctmp != "adsig")  stop("object is not of class adsig")
  if(is.null(ind)) stop("must provide name of individual to plot")
  if(sum(match(x$individuals,ind), na.rm=TRUE) != 1) stop("individual not recognised")
  if(sum(!is.na(match(x$individuals,popA))) != length(popA)) stop("individual in PopA not recognised")
 if(sum(!is.na(match(x$individuals,popB))) != length(popB)) stop("individual in PopB not recognised")
  
 #set up generic labels
  if(is.null(ylab)) ylab <- "signal"

 if(is.null(main)){ main <- paste("ID:", x$individuals[ind], ", Window size:" ,x$window.size,sep=" ")}


      #Set plot ranges
      tempP <- t(x$signals) #take signals
	#set x locations for plotting (genetic distance or SNP number)
	if(is.null(x$gendist)){
      loc <- 1:dim(tempP)[2]
	  if(is.null(xlab)) xlab <- "SNP number"}else{
	loc <- x$gendist
	  if(is.null(xlab)) xlab <- "genetic distance"
	}
      #range of ancestral populations
      if(!is.null(popA)){
        Amin <- apply(tempP[popA,],c(2),FUN=min)
        Amax <- apply(tempP[popA,],c(2),FUN=max)  }
      if(!is.null(popB)){
        Bmin <- apply(tempP[popB,],c(2),FUN=min)
        Bmax <- apply(tempP[popB,],c(2),FUN=max)}
      #set ylim
      if(is.null(ylim)){ ylim <- range(tempP, na.rm=TRUE)}

      
      #do the plot
      plot(loc,tempP[ind,], xlab=xlab, ylab=ylab,main=main, ylim=ylim,type="l",cex.main=1);abline(h=c(-1,0,1),lty=3,col=3);
      abline(h=mean(tempP[ind,],),col=2)
      if(!is.null(popB)){
        polygon(c(loc,rev(loc)), c((Bmin), rev(Bmax)), col = "lightblue", border = NA)}
      if(!is.null(popA)){
        polygon(c(loc,rev(loc)), c((Amax), rev(Amin)), col = "lightgreen", border = NA)}
      lines(loc,tempP[ind,],lwd=1,col=1)   
    		}
