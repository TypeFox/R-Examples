#
# vim:set ff=unix expandtab ts=2 sw=2:
deSolve.lsoda.wrapper=function(
### The function serves as a wrapper for lsoda using a much simpler interface which allows the use 
### of matrices in the definition of the derivative. 
### To use lsoda we have to convert our vectors to lists, define tolerances and so on.
### This function does this for us , so we don't need to bother about it.
	       t,	##<< A row vector containing the points in time where the solution is sought.
	       ydot,    ##<< The function of y and t that computes the derivative for a given 
	       ## point in time and a column vector y.
	       startValues ##<< A column vector with the starting values.
	       ){
   
   parms=NULL
   #my.atol <- 1e-6
   #rtol=1e-4
   lsexamp <- function(t, y,parms)
     {
	yv=cbind(y)
	YD=ydot(y,t)
	yd=as.vector(YD)
       #list(yd,c(massbalance=sum(y))) we could add other output parameter if we were interested
       list(yd)
     }
   #out <- lsoda(startValues,t,lsexamp, parms, rtol, atol= my.atol)
   out <- lsoda(startValues,t,lsexamp)
   #out <- rk(startValues,t,lsexamp,parms)
      #print(paste("out=",out))
      #print(out)
   # The output of lsoda is unsuiteable for our needs for two reasons
   # 1.) It also returns the time vector in column 1 
   # 2.) the columns get names instead of the default numbers created
   #     by the matrix function
   # we threrefore extract the information and store it in a new matrix witch will be t 
   n=length(startValues)
   if (n==1) { Yt=matrix(ncol=n,out[,-1])}
   else {Yt=out[,-1]}
   #print("Yt=")
   #print(Yt)
   #determine the number of pools 
   #determine the number of time values for which the solution is sought
   tn=length(t) 
   Y=matrix(ncol=n,nrow=length(t))
   #print(Yt[,1])
   for (i in 1:n){
      #print(paste("i=",i))
      Y[,i]=Yt[,i]
   }
   return(Y)
   ### A matrix. Every column represents a pool and every row a point in time
}
