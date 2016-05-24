## third.R
## estimate the third order moment
## Jan 2014
# inputs
#       data       - time-series
#       n.lag    - number of lags, max = length of time series
#       centre	  - centre series by subtracting mean (TRUE/FALSE)
#       outmax  - display lags of maxima and minima (TRUE/FALSE)
#       plot  - ggplot plot of the third order moment (TRUE/FALSE)


third=function(data,n.lag,centre=TRUE,outmax=TRUE,plot=TRUE){

  xaxis <- yaxis <- zaxis <- NULL # Setting some variables to NULL first (for R CMD check)
  
  nsamp=length(data);
  if(nsamp<10){cat('warning n<10\n')}
# ---------------- cumulants in non-redundant region -----------------
     XXX = matrix(data=0,n.lag+n.lag+1,n.lag+n.lag+1);
     if(centre==TRUE){centred=data-mean(data)}else{centred=data}
	 if(plot==TRUE){count=0}
     for (d in 0:n.lag){
         for (k in d:n.lag){
             large=max(c(d,k,0));
             XXX[d+n.lag+1,k+n.lag+1]=sum(centred[1:(nsamp-large)]*centred[(1+d):(nsamp-large+d)]*centred[(1+k):(nsamp-large+k)])/nsamp;
             XXX[n.lag+1+k,n.lag+1+d]=XXX[d+n.lag+1,k+n.lag+1]; # Symmetry
             XXX[n.lag+1-d,n.lag+1+k-d]=XXX[d+n.lag+1,k+n.lag+1]; # Symmetry
             XXX[n.lag+1+k-d,n.lag+1-d]=XXX[d+n.lag+1,k+n.lag+1]; # Symmetry
             XXX[n.lag+1+d-k,n.lag+1-k]=XXX[d+n.lag+1,k+n.lag+1]; # Symmetry
             XXX[n.lag+1-k,n.lag+1+d-k]=XXX[d+n.lag+1,k+n.lag+1]; # Symmetry
   	         if(plot==TRUE){
			    frame=data.frame(xaxis=d,yaxis=k,zaxis=XXX[d+n.lag+1,k+n.lag+1])
			    if(count==0){for.plot=frame}else{for.plot=rbind(for.plot,frame)}
				count=count+1
			 }
         }
     }
  
     waxis=-n.lag:n.lag;

# Lags of minima and maxima
  if(outmax==TRUE){
    cat('Maximum at (including symmetries)\n');
	cat(which(XXX==max(XXX),arr.ind=T)-n.lag-1,'\n')
    cat('Minimum at (including symmetries)\n');
	cat(which(XXX==min(XXX),arr.ind=T)-n.lag-1,'\n')
  }

# Lags of minima and maxima
  if(plot==TRUE){
     gplot=ggplot(for.plot, aes(xaxis, yaxis, z = zaxis))+
           stat_contour()+
           geom_tile(aes(fill = zaxis))
     print(gplot)
  }
  
  to.return=list()
  to.return$waxis=waxis
  to.return$third=XXX
  return(to.return)
 
} # end of function

