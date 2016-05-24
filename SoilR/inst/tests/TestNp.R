#
# vim:set ff=unix expandtab ts=2 sw=2:
TestNp=function(){
   require(RUnit)
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
   Yeuler	=TwopParallel(times,k1,k2,c10,c20,SoilR.euler)
   Yode	=TwopParallel(times,k1,k2,c10,c20,deSolve.lsoda.wrapper)
   checkEquals(target, Yode, "test numeric solution computed by deSolve.lsoda.wrapper function against analytical",
            #tolerance = .Machine$double.eps^0.5, 
            tolerance = tol,
   )
   checkEquals(target, Yeuler, "test numeric solution computed by the euler mehtod against analytical",
            tolerance = tol,
   	 )
   k3=-0.5
   c30=3
   target=ThreepParallelAnalytical(times,k1,k2,k3,c10,c20,c30)
   Yeuler=ThreepParallel(times,k1,k2,k3,c10,c20,c30,SoilR.euler)
   
   checkEquals(target, Yeuler, "test numeric solution computed by the euler mehtod against analytical",
            tolerance = tol,
   )
}
