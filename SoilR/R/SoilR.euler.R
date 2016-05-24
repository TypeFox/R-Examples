#
# vim:set ff=unix expandtab ts=2 sw=2:
SoilR.euler=function
### This function can solve arbitrary first order ode systems with the euler forward 
### method and an
### adaptive time-step size control given a tolerance for the deviation of a coarse and fine 
### estimate of the change in y for the next time step.
### It is an alternative to \code{\link{ode}} and has the same interface.
### It is much slower than ode and should probably be considered less capable in solving stiff ode systems.
### However it has one main advantage, which consists in its simplicity.
### It is quite easy to see what is going on inside it.
### Whenever you don't trust your implementation of another (more efficient but probably also more complex)
### ode solver, just compare the result to what this method computes.
(times,		##<< A row vector containing the points in time where the solution is sought.
 ydot,		##<< The function of y and t that computes the derivative for a given point in time and a column vector y.
 startValues		##<< A column vector with the initial values.
 ){
   inc=1.5
   tol=10**(-10) # this is the tolerance that is allowed for the difference of estimation and solution for the
   #next step
   ttol=10**(-10) #this is only a technical value that avoids an infinite loop caused by roundoferrors for time
   minstep=10**(-5) #the smalles timestep allowed to avoid locking

   #determine the number of unknowns
   n=nrow(startValues)
   #determine the number of time values for which the solution is sought
   tn=length(times) 
   #we give a startvalue for the timestepsize
   Y=matrix(nrow=n,ncol=tn)
   y=startValues
   Y[,1]=y
   #now devide the time in to steps 
      for (j in 2:tn){
	 targettime=times[j]
	 t=times[j-1]
	 stepsize=(targettime-t)/10
	 y=Y[,j-1] #we always start at the boundaries of timesteps
	 while (t< targettime-ttol){
	    #store in case we have to revert the next step
	    y0=y
	    dy=stepsize*ydot(y,t)
	    #make a prediction yp=y(t+stepsize) 
	    yp=y+dy
	    #now devide the stepsize in k subintervals and compute the value 
	    # fpr y(t+stepsize) in k steps
	    k=4
	    smallstep=stepsize/k
	    for (i in 1:k){
	       t=t+smallstep
	       dy=smallstep*ydot(y,t)
	       y=y+dy
	    }
	    #compare the results
	    ydiff=sum((yp-y)*(yp-y))
	    if (ydiff>tol){
	       # do it again
	       t=t-stepsize
	       y=y0
	       #in smaller steps
	       stepsize=stepsize/2
	    }
	    else{
	       # become bolder but look ahead a bit 
	       rest=targettime-t
	       planedstep=inc*stepsize
	       # if the target time is far enough ahead for the planed timestep
	       if (rest>planedstep)
		  #increase timestep but with a look at the target
		  #because we want to avoid a very small last timestep
		  if (rest-planedstep<minstep)
		     #since the planed step would bring us too near the target but not
		     #quite there we decide to suggest two steps instead with the
		     {stepsize=rest/2}
		  else
		     # or prooced as normal
		     {stepsize=stepsize*inc}
	       #otherwise make a direct leap to the next targettime
	       else{stepsize=rest}
	    }
	 }
	 Y[,j]=y
      }
   Yt=t(Y)
   return(Yt)
}


