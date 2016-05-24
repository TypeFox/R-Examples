t_start=1960
t_end=2010
tn=220
timestep=(t_end-t_start)/tn 
t=seq(t_start,t_end,timestep) 
n=3
At=new(Class="BoundLinDecompOp",
  t_start,
  t_end,
  function(t0){
        matrix(nrow=n,ncol=n,byrow=TRUE,
          c(-1,    0.1,    0, 
             0.5  , -0.4,    0,   
             0,    0.2,   -0.1)
        )
  }
) 
 
c0=c(100, 100, 100)
F0=ConstFc(c(0,10,10),"Delta14C")
#constant inputrate
inputFluxes=new(
  "TimeMap",
  t_start,
  t_end,
  function(t0){matrix(nrow=n,ncol=1,c(10,10,10))}
) 
# we have a dataframe representing the C_14 fraction 
# note that the time unit is in years and the fraction is given in
# the Absolute Fraction Modern format.
# This means that all the other data provided are assumed to have the same value
# This is especially true for the decay constants to be specified later
inputFc=BoundFc(C14Atm_NH,format="Delta14C")
# add the C14 decay to the matrix which is done by a diagonal matrix which does not vary over time
# we assume a half life th=5730 years
th=5730
k=log(0.5)/th #note that k is negative and has the unit y^-1

mod=GeneralModel_14(t=t,A=At,ivList=c0,initialValF=F0,inputFluxes=inputFluxes,inputFc=inputFc,di=k)
#start plots
par(mfrow=c(3,2))
   lt1=1;  lt2=2; lt3=3 
   col1=1;  col2=2; col3=3
   # plot the C and C14 curves
   Ct=getC(mod)
   plot(t,Ct[,1],type="l",lty=lt1,col=col1, ylim=c(0,200),
        ylab="C stocks (arbitrary units)",xlab="Time") 
   lines(t,Ct[,2],type="l",lty=lt2,col=col2) 
   lines(t,Ct[,3],type="l",lty=lt3,col=col3,lwd=2) 
   legend(
      "topright",
      c("C in pool 1",
        "C in pool 2",
        "C in pool 3"
      ),
      lty=c(lt1,lt2,lt3),
      col=c(col1,col2,col3)
   )
   C14t=getC14(mod)
   plot(t,C14t[,1]/1000,type="l",lty=lt1,col=col1, ylim=c(0,100),
        ylab="14C stocks (arbitrary units)",xlab="Time") 
   lines(t,C14t[,2]/1000,type="l",lty=lt2,col=col2) 
   lines(t,C14t[,3]/1000,type="l",lty=lt3,col=col3) 
   legend(
      "topright",
      c("14C in pool 1",
        "14C in pool 2",
        "14C in pool 3"
      ),
      lty=c(lt1,lt2,lt3),
      col=c(col1,col2,col3)
   )
   #now plot the C14 Fraction in the atmosphere and compute the C14/C fraction of in the Soil 

   FC14=getF14(mod)
   plot(C14Atm_NH, type="l",xlim=c(1960,2010))
   lines(t,FC14[,1],lty=lt1,col=col1) 
   lines(t,FC14[,2],lt2,type="l",lty=lt2,col=col2) 
   lines(t,FC14[,3],type="l",lty=lt3,col=col3) 
   legend("topleft",c(
                      expression(F[1]),
                      expression(F[2]),
                      expression(F[3])
                    )
   ,lty=c(lt1,lt2,lt3),col=c(col1,col2,col3))


#now compute the release flux
   Rt=getReleaseFlux(mod)
   plot(
     t,
     Rt[,1],
     type="l",
     lty=lt1,
     col=col1,
     ylab="C Release Flux (arbitrary units)",
     xlab="Time",
     ylim=c(0,50)
   ) 
   lines(t,Rt[,2],lt2,type="l",lty=lt2,col=col2) 
   lines(t,Rt[,3],type="l",lty=lt3,col=col3) 
   legend("topleft",c("RF1","RF2","RF3"),lty=c(lt1,lt2,lt3),col=c(col1,col2,col3))
   #now compute the c14 release flux
   R14t=getReleaseFlux14(mod)/1000
   plot(
     t,
     R14t[,1],
     type="l",
     lty=lt1,
     col=col1,
     ylab="C14 Release Flux (arbitrary units)",
     xlab="Time"
   ) 
   lines(t,R14t[,2],lt2,type="l",lty=lt2,col=col2) 
   lines(t,R14t[,3],type="l",lty=lt3,col=col3) 
   legend("topleft",c(
                      expression(RF[14]^1),
                      expression(RF[14]^2),
                      expression(RF[14]^3)
                    )
   ,lty=c(lt1,lt2,lt3),col=c(col1,col2,col3))

R14m=getF14R(mod)
C14m=getF14C(mod)
plot(C14Atm_NH, type="l",xlim=c(1960,2010),col=4)
lines(t,C14m) 
lines(t,R14m,col=2) 
legend(
  "topright",
  c("Atmosphere","Mean SOM-14C","Mean Release 14C"),
  lty=rep(1,3),
  col=c(4,1,2),
  bty="n"
)
par(mfrow=c(1,1))
