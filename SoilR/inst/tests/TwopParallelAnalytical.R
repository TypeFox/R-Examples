#
# vim:set ff=unix expandtab ts=2 sw=2:
TwopParallelAnalytical<-structure(function
### The function computes the solution for two independent (parallel) pools analytically
(
 ### title<< Simple function arguments
 time,	##<< a vector containing the points in time where the solution is sougth
 k1,    ##<< the rate of decay in pool 1
 k2,	##<< the same for pool 2
 c10,	##<< initial concentration for pool 1
 c20 	##<< the same for pool 2
 ){
   ##details<< 
   ##Each row of the matrix is computed independently
   ##describe<<
      y1=c10*exp(k1*time) ##<<values for pool one
      y2=c20*exp(k2*time) ##<<values for pool two
      Y=matrix(ncol=2,nrow=length(time))
      Y[,1]=y1
      Y[,2]=y2
   return( Y )
   ##note<< An executable example is attached to the function as attribute "ex"
   ## you can extract it with attr(TwopParrallelAnalytical,"ex")

   ##seealso<< \code{\link{TwopParallelAnalytical}}
}
,
#ex=TestNp
ex=function(){
   t_start=0 
   t_end=10 
   tn=10 
   tol=.02/tn
   print(tol)
   timestep=(t_end-t_start)/tn 
   times=seq(t_start,t_end,timestep) 
   k1=-0.5 
   k2=-0.2
   c10=1
   c20=2
   target=TwopParallelAnalytical(times,k1,k2,c10,c20)
   eulerMod=TwopParallelModel(times,k1,k2,c10,c20,SoilR.euler)
   Yeuler=getCcontent(eulerMod)
   odeMod=TwopParallelModel(times,k1,k2,c10,c20,deSolve.lsoda.wrapper)
   Yode=getCcontent(odeMod)
   checkEquals(target, Yode, "test numeric solution computed by deSolve.lsoda.wrapper function against analytical",
            #tolerance = .Machine$double.eps^0.5, 
            tolerance = tol,
   )
   checkEquals(target, Yeuler, "test numeric solution computed by the euler mehtod against analytical",
            tolerance = tol,
   	 )
}
)


