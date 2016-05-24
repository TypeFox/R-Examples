# This test function is automatically produced by the python script:/home/mm/SoilR/RPackages/SoilR/pkg/inst/tests/automatic/Rexample.py
test.TwopParallel_ZeroDecayInputOnly=function(){
   require(RUnit)
   t_start=0
   t_end=2
   tn=100
   tol=.02/tn
   print(tol)
   timestep=(t_end-t_start)/tn
   t=seq(t_start,t_end,timestep)
   A=new("ConstLinDecompOp",matrix(
     nrow=2,
     ncol=2,
     c(
        0,  0,  
        0,  0
     )
   ))
   c01=3
   c02=2
   inputrates=new("TimeMap",t_start,t_end,function(t){return(matrix(
     nrow=2,
     ncol=1,
     c(
        0.100000000000000,  0.200000000000000
     )
   ))})
   Y=matrix(ncol=2,nrow=length(t))
   Y[,1]=c01 + 0.1*t
   Y[,2]=c02 + 0.2*t
   R=matrix(ncol=2,nrow=length(t))
   R[,1]=0
   R[,2]=0
   mod=GeneralModel(
    t,
    A,
    c(
       c01,
       c02
    ),
   inputrates,
   deSolve.lsoda.wrapper
   )
   Yode=getC(mod) 
   Rode=getReleaseFlux(mod) 
#begin plots 
   lt1=2
   lt2=4
   pdf(file="runit.TwopParallel_ZeroDecayInputOnly.pdf",paper="a4")
   m=matrix(c(1,2,3,4),4,1,byrow=TRUE)
   layout(m)
   plot(t,Y[,1],type="l",lty=lt1,col=1,ylab="Concentrations",xlab="Time")
   lines(t,Yode[,1],type="l",lty=lt2,col=1)
   lines(t,Y[,2],type="l",lty=lt1,col=2)
   lines(t,Yode[,2],type="l",lty=lt2,col=2)
   legend(
   "topright",
     c(
     "anlytic sol for pool 1",
     "numeric sol for pool 1",
     "anylytic sol for pool 2",
     "numeric sol for pool 2"
     ),
     lty=c(lt1,lt2),
     col=c(1,1,2,2)
   )
   plot(t,R[,1],type="l",lty=lt1,col=1,ylab="Respirationfluxes",xlab="Time",ylim=c(min(R),max(R)))
   lines(t,Rode[,1],type="l",lty=lt2,col=1)
   lines(t,R[,2],type="l",lty=lt1,col=2)
   lines(t,Rode[,2],type="l",lty=lt2,col=2)
   legend(
   "topright",
     c(
     "anlytic sol for pool 1",
     "numeric sol for pool 1",
     "anylytic sol for pool 2",
     "numeric sol for pool 2"
     ),
     lty=c(lt1,lt2),
     col=c(1,1,2,2)
   )
   dev.off()
# end plots 
# begin checks 
   checkEquals(
    Y,
    Yode,
    "test numeric solution for C-Content computed by the ode mehtod against analytical",
    tolerance = tol,
   )
   checkEquals(
    R,
    Rode,
    "test numeric solution for Respiration computed by the ode mehtod against analytical",
    tolerance = tol,
   )

 }