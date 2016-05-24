# This test function is automatically produced by the python script:/home/mm/SoilR/RPackages/SoilR/pkg/inst/tests/automatic/Rexample.py
test.op=function(){
   require(RUnit)
   c1=1
   k1=1/10
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
        -k1
     )
   ))
   inputrates=new("TimeMap",t_start,t_end,function(t){return(matrix(
     nrow=1,
     ncol=1,
     c(
        0
     )
   ))})
   Y=matrix(ncol=1,nrow=length(t))
   Y[,1]=c1*exp(-t/10)
   R=matrix(ncol=1,nrow=length(t))
   R[,1]=c1*exp(-t/10)/10
meanTransitTime=1/k1
   mod=GeneralModel(
    t,
    A,
    c(
       c1
    ),
   inputrates,
   deSolve.lsoda.wrapper
   )
   Yode=getC(mod) 
   Rode=getReleaseFlux(mod) 
   meanTransitTimeode=getMeanTransitTime(
        A,
    c(
       c1
    )
)
   TTDode=getTransitTimeDistributionDensity(
        A,
    c(
       c1
    )
,t
)
#begin plots 
   lt1=2
   lt2=4
   pdf(file="runit.op.pdf",paper="a4")
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
   plot(t,TTDode,type="l",lty=lt1,col=1,ylab="TransitTimeDistributionDensity",xlab="Time")
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
   checkEquals(
    meanTransitTime,
    meanTransitTimeode,
    "test numeric solution for the mean transit Tiye computed by the ode mehtod against analytical value taken from manzoni et al",
    tolerance = tol,
   )

 }