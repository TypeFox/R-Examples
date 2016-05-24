# This test function is automatically produced by the python script:/home/mm/SoilR/RPackages/SoilR/pkg/inst/tests/automatic/Rexample.py
test.TwopFeedback=function(){
   require(RUnit)
   c1=1
   c2=0
   r=1/4
   k1=1/10
   k2=1/5
   f=1
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
        -k1,  k1*(-r + 1),  
        f*k2,  -k2
     )
   ))
   inputrates=new("TimeMap",t_start,t_end,function(t){return(matrix(
     nrow=2,
     ncol=1,
     c(
        0,  0
     )
   ))})
   Y=matrix(ncol=2,nrow=length(t))
   Y[,1]=c1*(4*exp(t*(-3/20 - sqrt(7)/20))/((2/3 + 2*sqrt(7)/3)*(1 + sqrt(7))*(-sqrt(7)/3 + 7/3)) - 4*(-4/((2/3 + 2*sqrt(7)/3)**2*(1 + sqrt(7))*(-sqrt(7)/3 + 7/3)) + 1/(2/3 + 2*sqrt(7)/3))*exp(t*(-3/20 + sqrt(7)/20))/(-sqrt(7) + 1)) + c2*(-4*exp(t*(-3/20 - sqrt(7)/20))/((1 + sqrt(7))*(-sqrt(7)/3 + 7/3)) - 16*exp(t*(-3/20 + sqrt(7)/20))/((2/3 + 2*sqrt(7)/3)*(1 + sqrt(7))*(-sqrt(7) + 1)*(-sqrt(7)/3 + 7/3)))
   Y[,2]=c1*(-exp(t*(-3/20 - sqrt(7)/20))/((2/3 + 2*sqrt(7)/3)*(-sqrt(7)/3 + 7/3)) + (-4/((2/3 + 2*sqrt(7)/3)**2*(1 + sqrt(7))*(-sqrt(7)/3 + 7/3)) + 1/(2/3 + 2*sqrt(7)/3))*exp(t*(-3/20 + sqrt(7)/20))) + c2*(exp(t*(-3/20 - sqrt(7)/20))/(-sqrt(7)/3 + 7/3) + 4*exp(t*(-3/20 + sqrt(7)/20))/((2/3 + 2*sqrt(7)/3)*(1 + sqrt(7))*(-sqrt(7)/3 + 7/3)))
   R=matrix(ncol=2,nrow=length(t))
   R[,1]=c1*(4*exp(t*(-3/20 - sqrt(7)/20))/((2/3 + 2*sqrt(7)/3)*(1 + sqrt(7))*(-sqrt(7)/3 + 7/3)) - 4*(-4/((2/3 + 2*sqrt(7)/3)**2*(1 + sqrt(7))*(-sqrt(7)/3 + 7/3)) + 1/(2/3 + 2*sqrt(7)/3))*exp(t*(-3/20 + sqrt(7)/20))/(-sqrt(7) + 1))/40 + c2*(-4*exp(t*(-3/20 - sqrt(7)/20))/((1 + sqrt(7))*(-sqrt(7)/3 + 7/3)) - 16*exp(t*(-3/20 + sqrt(7)/20))/((2/3 + 2*sqrt(7)/3)*(1 + sqrt(7))*(-sqrt(7) + 1)*(-sqrt(7)/3 + 7/3)))/40
   R[,2]=0
meanTransitTime=(k1*(-r + 1) + k2)/(k1*k2*r)
   mod=GeneralModel(
    t,
    A,
    c(
       c1,
       c2
    ),
   inputrates,
   deSolve.lsoda.wrapper
   )
   Yode=getC(mod) 
   Rode=getReleaseFlux(mod) 
   meanTransitTimeode=getMeanTransitTime(
        A,
    c(
       c1,
       c2
    )
)
   TTDode=getTransitTimeDistributionDensity(
        A,
    c(
       c1,
       c2
    )
,t
)
#begin plots 
   lt1=2
   lt2=4
   pdf(file="runit.TwopFeedback.pdf",paper="a4")
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