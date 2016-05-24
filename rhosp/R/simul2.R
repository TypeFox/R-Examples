"simul2" <-
function(nbPatient,disXi,disYi,toplot=TRUE)
{
#disXi is a 3 elements list : rangen <=>  a random positive variable generator
#       ; nbparam <=> number of parameter of this distribution and param <=> a list of parameters
#disYi is a 3 elements list : rangen <=> a random positive variable generator
#       ; nbparam <=> number of parameter of this distribution and param <=> a list of parameters

       	#simulation of variables Xi, Yi and Zi
       	if(disXi$nbparam==1) x<-disXi$rangen(nbPatient,disXi$param[[1]]);
       	if(disXi$nbparam==2) x<-disXi$rangen(nbPatient,disXi$param[[1]],disXi$param[[2]]);
       	if(disXi$nbparam==3) x<-disXi$rangen(nbPatient,disXi$param[[1]],disXi$param[[2]],disXi$param[[3]]);
       
       	if(disYi$nbparam==1) y<-disYi$rangen(nbPatient,disYi$param[[1]]);
       	if(disYi$nbparam==2) y<-disYi$rangen(nbPatient,disYi$param[[1]],disYi$param[[2]]);
       	if(disYi$nbparam==3) y<-disYi$rangen(nbPatient,disYi$param[[1]],disYi$param[[2]],disYi$param[[3]]);
       
	   	xCumsum<-cumsum(x);
	   	yCumsum<-cumsum(y);
			
	 
	   	z<-array(0,nbPatient);
		T<-array(0,1);
		
		for (i in 1:nbPatient)
		{
			if(i>1)
			{
				if(y[i]<x[i])
				{
					z[i]<-1;
				}
				else
				{
					z[i]<-0;
				}
			}
		}
	
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

