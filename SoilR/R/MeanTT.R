#
# vim:set ff=unix expandtab ts=2 sw=2:
MeanTT=function# mean transit time for a general one pool model
### The function computes the mean transit time for one pool of a possibly nonlinear model
(OdotLin,	##<< The outputrate of this pool as a linear operator (a function of y and t)
times		##<< A vector containing the points in time where the solution is sought.
){

   	 tstart=min(times)
   	 tend=max(times)
   	 ######################################################################
   	 # We assume  the amount of matter in the pool
   	 # at tstart=0 C0=1  ignore all inputs and compute the time dependent amount 
   	 # still in the system after time t C(t)
   	 # the probality for a particle to have a traveling time 
   	 # through the system of at least t 
   	 # is the given by the amount left C(t)/C0 =C(t) (since C0=1)
   	 impulsiveInputSolution=splinefun(times,solver(times,OdotLin,1))

	 c=c("black","red","green","blue")
   	 lts=c(1,2)
   	 lws=c(8,4)
   	 ## The impulsive InputSolution is what we would have observed when we had dyed (or marked)
   	 ## part of the Input at time 0
   	 ## We can treat it as  representative for any amount of input however tiny 
   	 ## since we have a linear system.
   	 ## For t in [ 0 ,tend]  we now know the probability of having traveled longer
   	 ## than t
   	 ## p=int_T^inf phi(t) dt
   	 ## phi(T)= -d/dT p=-OdotLin(Y(t),t)
   	 phi=function(t){-OdotLin(impulsiveInputSolution(t),t)}
   	 phis=splinefun(times,sapply(times,phi))
	 plot(times,phis(times))
   	 #lines(times,phi(times),col=c[4],lty=lts[1])#,ylim=c(0,1))
   	 #pf("phi(tstart)")
   	 phiOde=function(Y,t){phis(t)}
   	 p=function(ta){integrate(phis,lower=ta,upper=tend)$value}
   	 plot(times,impulsiveInputSolution(times),type="l",lwd=lws[1],col=c[1],lty=lts[1])
   	 lines(times,sapply(times,p),lwd=lws[2],col=c[2],lty=lts[2],ylim=c(0,1))
   	 plot(times,sapply(times,p),lwd=lws[2],col=c[2],lty=lts[2],ylim=c(0,1))
	 legend(
   	 "topright",
   	   c( "solution of linear system for impulsive Input","P(t<age<tmax)"),
   	   lty=lts,
   	   col=c(c[1],c[2])
   	 )
   	 #lines(times,impulsiveInputSolution(times),col=c[2],lty=lts[2])
   	 ## to get the mean transit time we have to multiply the density with the age
   	 ######################################################################
   	 phiTT=function(Y,tt){phiOde(Y,tt)*tt}
   	 mean_TT=solver(times,phiTT,0)
   	 plot(times,mean_TT,col=c[4],lty=lts[1])#,ylim=c(0,1))
	 return(mean_TT)
   ### A vector containing the mean transit time for the specified times
   }
