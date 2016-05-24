"simul" <-
function(nbPatient,disXi,disP,toplot=TRUE)
{
#disXi is a 3 elements list : rangen <=>  a random positive variable generator
#       ; nbparam <=> number of parameter of this distribution and param <=> a list of parameters
#disP is a 3 elements list : disfun <=> a distribution function
#       ; nbparam <=> number of parameter of this distribution and param <=> a list of parameters

       #simulation of variables Xi and Zi
       if(disXi$nbparam==1) x<-disXi$rangen(nbPatient,disXi$param[[1]]);
       if(disXi$nbparam==2) x<-disXi$rangen(nbPatient,disXi$param[[1]],disXi$param[[2]]);
       if(disXi$nbparam==3) x<-disXi$rangen(nbPatient,disXi$param[[1]],disXi$param[[2]],disXi$param[[3]]);
	   xCumsum<-cumsum(x);
	
	   p_x<-0;
	
	   if(disP$nbparam==1) p_x<-sapply(x,disP$param[[1]],FUN=disP$disfun);
       if(disP$nbparam==2) p_x<-sapply(x,disP$param[[1]],disP$param[[2]],FUN=disP$disfun);
       if(disP$nbparam==3) p_x<-sapply(x,disP$param[[1]],disP$param[[2]],disP$param[[3]],FUN=disP$disfun);
       if(disP$nbparam==4) p_x<-sapply(x,disP$param[[1]],disP$param[[2]],disP$param[[3]],disP$param[[4]],FUN=disP$disfun);

	   z<-rbinom(nbPatient,1,p_x);
	   T<-array(0,1);
	
	   #plot
	   if(toplot) plot(xCumsum,z,type="h",xlab="time",ylab="side effect");
	
	   #calculus of T
	   iPremVidT<-1;
	   iAncienT<-1;
	   for (i in 1:nbPatient)
	   {
	           if(z[i]==1)
	           {
	                   if(iPremVidT==1)
	                   {
	                           T[iPremVidT]<-xCumsum[i];
	                           iPremVidT<-iPremVidT+1;
	                           iAncienT<-i;
	                   }
	                   else
	                   {
	                           T[iPremVidT]<-xCumsum[i]-xCumsum[iAncienT];
	                           iPremVidT<-iPremVidT+1;
	                           iAncienT<-i;
	                   }
	           }
	   }
	   return (T);
}

