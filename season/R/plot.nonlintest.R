### Plot the results of the non-linear test

plot.nonlintest<-function(x, ...){

  z <- NULL # Setting some variables to NULL first (for R CMD check)
    
  ## Check
  if (class(x) != "nonlintest"){
    stop("Object must be of class 'nonlintest'")
  } 
  
if(max(abs(x$region))==0){cat('No points of the third-order moment exceed the test limits\n')}

  # Plot (only points that exceed limits)
if(max(abs(x$region))!=0){
   counter=0;
   for(r in 0:x$n.lag){
      for(s in r:x$n.lag){
         frame=data.frame(r=r,s=s,z=x$region[r+1,s+1]);
         if(counter==0){for.plot=frame}else{for.plot=rbind(for.plot,frame)}
		 counter=counter+1;
      }
   }

     gplot=ggplot(for.plot, aes(r, s, z = z))+
           stat_contour()+
           geom_tile(aes(fill = z))
     print(gplot)
 #zlab='Area outside the test limits', ...)
 #  mtext(side=3,'Points of 3rd order moment that exceed limits');
} # end of if
} # end of function

