# This test function is automatically produced by the python script:/home/mm/SoilR/RPackages/SoilR/pkg/inst/tests/automatic/Rexample.py
test.OnePool_C14_ZeroDecay_constantInputrate_c14_fromDelta14C=function(){
   require(RUnit)
   t_start=0
   t_end=2
   tn=100
   tol=.02/tn
   print(tol)
   timestep=(t_end-t_start)/tn
   t=seq(t_start,t_end,timestep)
   A=new("ConstLinDecompOp",matrix(
     nrow=1,
     ncol=1,
     c(
        0
     )
   ))
   c01=1
   inputrates=new("TimeMap",t_start,t_end,function(t){return(matrix(
     nrow=1,
     ncol=1,
     c(
        1
     )
   ))})
   f01=1
   initialF=ConstFc(    c(
       f01
    ),
 format="Delta14C")
   Fc=BoundFc(function(t){0.5},t_start,t_end,format="Delta14C")
   th=5730
   k=log(0.5)/th
   Y=matrix(ncol=1,nrow=length(t))
   Y[,1]=c01 + t
   R=matrix(ncol=1,nrow=length(t))
   R[,1]=0
   Y14=matrix(ncol=1,nrow=length(t))
   Y14[,1]=c01*(0.001*f01 + 1.0)*exp(-t*log(2)/5730) + 5732.865/log(2) - 5732.865*exp(-t*log(2)/5730)/log(2)
   F14=matrix(ncol=1,nrow=length(t))
   F14[,1]=-1000 + 1000*(c01*(0.001*f01 + 1.0)*exp(-t*log(2)/5730) + 5732.865/log(2) - 5732.865*exp(-t*log(2)/5730)/log(2))/(c01 + t)
   mod=GeneralModel_14(
    t=t,
    A=A,
ivList=    c(
       c01
    ),
initialValF=   initialF,
inputFluxes=   inputrates,
inputFc=   Fc,
di=   k,
solverfunc=   deSolve.lsoda.wrapper
   )
   Y14ode=getC14(mod) 
   F14ode=getF14(mod) 
   Yode=getC(mod) 
   Rode=getReleaseFlux(mod) 
#begin plots 
   lt1=2
   lt2=4
   pdf(file="runit.OnePool_C14_ZeroDecay_constantInputrate_c14_fromDelta14C.pdf",paper="a4")
   m=matrix(c(1,2,3,4),4,1,byrow=TRUE)
   layout(m)
   plot(t,Y[,1],type="l",lty=lt1,col=1,ylab="Concentrations",xlab="Time")
   lines(t,Yode[,1],type="l",lty=lt2,col=1)
   legend(
   "topright",
     c(
     "anylytic sol for pool 1",
     "numeric sol for pool 1"
     ),
     lty=c(lt1,lt2),
     col=c(1,1)
   )
   plot(t,R[,1],type="l",lty=lt1,col=1,ylab="Respirationfluxes",xlab="Time",ylim=c(min(R),max(R)))
   lines(t,Rode[,1],type="l",lty=lt2,col=1)
   legend(
   "topright",
     c(
     "anylytic sol for pool 1",
     "numeric sol for pool 1"
     ),
     lty=c(lt1,lt2),
     col=c(1,1)
   )
   plot(t,Y14[,1],type="l",lty=lt1,col=1,ylab="14C-Concentrations",xlab="Time",ylim=c(min(Y14),max(Y14)))
   lines(t,Y14ode[,1],type="l",lty=lt2,col=1)
   plot(t,F14[,1],type="l",lty=lt1,col=1,ylab="14C-C ratio ",xlab="Time",ylim=c(min(F14,F14ode),max(F14,F14ode)))
   lines(t,F14ode[,1],type="l",lty=lt2,col=1)
   legend(
   "topright",
     c(
     "anylytic sol for pool 1",
     "numeric sol for pool 1"
     ),
     lty=c(lt1,lt2),
     col=c(1,1)
   )
   dev.off()
# end plots 
# begin checks 
   tol=.02*max(Y14)/tn
   checkEquals(
    Y14,
    Y14ode,
    "test numeric solution for 14C-Content computed by the ode mehtod against analytical",
    tolerance = tol,
   )
   checkEquals(
    F14,
    F14ode,
    "test numeric solution for F14 computed by the ode mehtod against analytical",
    tolerance = tol,
   )

 }