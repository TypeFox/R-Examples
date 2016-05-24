#
# vim:set ff=unix expandtab ts=2 sw=2:
ThreepParallelAnalytical=structure( function
### The function computes the solution for two independent (parallel) pools analytically
(
 ### title<< Simple function arguments
 time,	##<< a vector containing the points in time where the solution is sougth
 k1,    ##<< the rate of decay in pool 1
 k2,	##<< the same for pool 2
 k3,	##<< the same for pool 3
 c10,	##<< initial concentration for pool 1
 c20, 	##<< the same for pool 2
 c30 	##<< the same for pool 3
 ){
   ##details<< 
   ##Each row of the matrix is computed independently
   ##describe<<
      Y=matrix(ncol=3,nrow=length(time))
      Y[,1]=c10*exp(k1*time) ##<<values for pool one
      Y[,2]=c20*exp(k2*time) ##<<values for pool two
      Y[,3]=c30*exp(k3*time) ##<<values for pool two
   return( Y )
##seealso<< \code{\link{TwopParrallelAnalytical}}
##seealso<< \code{\link{ThreepParrallelModel}}
}

#,ex= TestNp
,ex=function(){
   t_start=0 
   t_end=10 
   tn=10 
   tol=.02/tn
   print(tol)
   timestep=(t_end-t_start)/tn 
   times=seq(t_start,t_end,timestep) 
   k1=-0.5 
   k2=-0.2
   k3=-0.2
   c10=1
   c20=2
   c30=3
   target=ThreepParallelAnalytical(times,k1,k2,k3,c10,c20,c30)
   
   eulerMod=ThreepParallelModel(times,k1,k2,k3,c10,c20,c30,SoilR.euler)
   Yeuler=getCcontent(eulerMod)
   
   odeMod=ThreepParallelModel(times,k1,k2,k3,c10,c20,c30,deSolve.lsoda.wrapper)
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


