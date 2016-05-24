# This test function is automatically produced by the python script:/home/mm/SoilR/RPackages/SoilR/pkg/inst/tests/automatic/Rexample.py
test.ThreepSerial_1=function(){
   require(RUnit)
   t_start=0
   t_end=2
   tn=100
   tol=.02/tn
   print(tol)
   timestep=(t_end-t_start)/tn
   t=seq(t_start,t_end,timestep)
   A=new("ConstLinDecompOp",matrix(
     nrow=3,
     ncol=3,
     c(
        -1/2,  1/2,  0,  
        0,  -1/3,  1/6,  
        0,  0,  -1
     )
   ))
   c01=3
   c02=2
   c03=1
   inputrates=new("TimeMap",t_start,t_end,function(t){return(matrix(
     nrow=3,
     ncol=1,
     c(
        1,  2,  3
     )
   ))})
   Y=matrix(ncol=3,nrow=length(t))
   Y[,1]=c01*exp(-t/2) + 2 - 2*exp(-t/2)
   Y[,2]=c01*(-3*exp(-t/2) + 3*exp(-t/3)) + c02*exp(-t/3) + 9 + 6*exp(-t/2) - 15*exp(-t/3)
   Y[,3]=c01*(exp(-t)/4 - exp(-t/2) + 3*exp(-t/3)/4) + c02*(-exp(-t)/4 + exp(-t/3)/4) + c03*exp(-t) + 9/2 - 11*exp(-t)/4 + 2*exp(-t/2) - 15*exp(-t/3)/4
   R=matrix(ncol=3,nrow=length(t))
   R[,1]=0
   R[,2]=c01*(-3*exp(-t/2) + 3*exp(-t/3))/6 + c02*exp(-t/3)/6 + 3/2 + exp(-t/2) - 5*exp(-t/3)/2
   R[,3]=c01*(exp(-t)/4 - exp(-t/2) + 3*exp(-t/3)/4) + c02*(-exp(-t)/4 + exp(-t/3)/4) + c03*exp(-t) + 9/2 - 11*exp(-t)/4 + 2*exp(-t/2) - 15*exp(-t/3)/4
   mod=GeneralModel(
    t,
    A,
    c(
       c01,
       c02,
       c03
    ),
   inputrates,
   deSolve.lsoda.wrapper
   )
   Yode=getC(mod) 
   Rode=getReleaseFlux(mod) 
#begin plots 
   lt1=2
   lt2=4
   pdf(file="runit.ThreepSerial_1.pdf",paper="a4")
   m=matrix(c(1,2,3,4),4,1,byrow=TRUE)
   layout(m)
   plot(t,Y[,1],type="l",lty=lt1,col=1,ylab="Concentrations",xlab="Time")
   lines(t,Yode[,1],type="l",lty=lt2,col=1)
   lines(t,Y[,2],type="l",lty=lt1,col=2)
   lines(t,Yode[,2],type="l",lty=lt2,col=2)
   lines(t,Y[,3],type="l",lty=lt1,col=3)
   lines(t,Yode[,3],type="l",lty=lt2,col=3)
   legend(
   "topright",
     c(
     "anlytic sol for pool 1",
     "numeric sol for pool 1",
     "anlytic sol for pool 2",
     "numeric sol for pool 2",
     "anylytic sol for pool 3",
     "numeric sol for pool 3"
     ),
     lty=c(lt1,lt2),
     col=c(1,1,2,2,3,3)
   )
   plot(t,R[,1],type="l",lty=lt1,col=1,ylab="Respirationfluxes",xlab="Time",ylim=c(min(R),max(R)))
   lines(t,Rode[,1],type="l",lty=lt2,col=1)
   lines(t,R[,2],type="l",lty=lt1,col=2)
   lines(t,Rode[,2],type="l",lty=lt2,col=2)
   lines(t,R[,3],type="l",lty=lt1,col=3)
   lines(t,Rode[,3],type="l",lty=lt2,col=3)
   legend(
   "topright",
     c(
     "anlytic sol for pool 1",
     "numeric sol for pool 1",
     "anlytic sol for pool 2",
     "numeric sol for pool 2",
     "anylytic sol for pool 3",
     "numeric sol for pool 3"
     ),
     lty=c(lt1,lt2),
     col=c(1,1,2,2,3,3)
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