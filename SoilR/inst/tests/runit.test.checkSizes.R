#
# vim:set ff=unix expandtab ts=2 sw=2:
test.checkSizes_AunequalV=function(){
   t_start=0 
   t_end=10 
   tn=50
   timestep=(t_end-t_start)/tn 
   t=seq(t_start,t_end,timestep) 
   A=TimeMap.new(
    t_start,
    t_end,
    function(times){matrix(nrow=3,ncol=3,byrow=TRUE,
        c(-1,    0,    0, 
         0.5,   -2,    0,   
           0,    1, -0.5)
   )    
   }
  )  
  I=TimeMap.new(
     t_start,
     t_end,
     function(times){
       # the matrix dimension does not agree with the matrix defined
       # in A
       matrix(nrow=2,ncol=1,byrow=TRUE,
           c(-1,    0)
       )
     }
  )
  checkException(checkSizes(A,I), msg="checkSizes should have stopped executions, but has not")
}
