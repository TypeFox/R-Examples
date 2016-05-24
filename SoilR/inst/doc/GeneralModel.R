### R code from vignette source 'GeneralModel.Rnw'

###################################################
### code chunk number 1: GeneralModel.Rnw:106-107
###################################################
library("SoilR")


###################################################
### code chunk number 2: GeneralModel.Rnw:111-124
###################################################
n=3;
t_start=1;t_end=2
At=BoundLinDecompOp(
  function(t0){
        matrix(nrow=n,ncol=n,byrow=TRUE,
          c(-0.39,    0,    0, 
             0.1, -0.35,    0,   
             0,    1/3,   -0.33)
        )
  },
  t_start,
  t_end
) 


###################################################
### code chunk number 3: GeneralModel.Rnw:130-135
###################################################
inputFluxes=BoundInFlux(
    function(t0){matrix(nrow=3,ncol=1,c(0.05,0,0))},
    t_start,
    t_end
)


###################################################
### code chunk number 4: GeneralModel.Rnw:139-142
###################################################
tn=500
timestep=(t_end-t_start)/tn
t=seq(t_start,t_end,timestep) 


###################################################
### code chunk number 5: GeneralModel.Rnw:147-148
###################################################
c0=c(0.5, 0.5, 0.5) 


###################################################
### code chunk number 6: GeneralModel.Rnw:151-152
###################################################
mod=GeneralModel(t,At,c0,inputFluxes) 


###################################################
### code chunk number 7: GeneralModel.Rnw:156-157
###################################################
Y_c=getC(mod) 


###################################################
### code chunk number 8: GeneralModel.Rnw:165-185
###################################################
lt1=1;  lt2=2; lt3=3 
col1=1; col2=2; col3=3 
plot(t,
	Y_c[, 1],
	type="l",
	lty=lt1,
	col=col1,
	ylab="C stocks (arbitrary units)",
	xlab="Time",
	ylim=c(min(Y_c),max(Y_c))
)
lines(t,Y_c[,2],type="l",lty=lt2,col=col2)
lines(t,Y_c[,3],type="l",lty=lt3,col=col3) 
legend(
    "topright",
    c("C in pool 1", "C in pool 2", "C in pool 3"),
    lty=c(lt1,lt2,lt3),
    col=c(col1,col2,col3),
    bty="n"
)


###################################################
### code chunk number 9: GeneralModel.Rnw:196-197
###################################################
Y_rf=getReleaseFlux(mod)


###################################################
### code chunk number 10: GeneralModel.Rnw:201-208
###################################################
    plot(t,Y_rf[,1],type="l",lty=lt1,col=col1,
    ylab="C Release (arbitrary units)",
    xlab="Time", ylim=c(0,0.2))
    lines(t,Y_rf[,2],lt2,type="l",lty=lt2,col=col2)
    lines(t,Y_rf[,3],type="l",lty=lt3,col=col3)
    legend("topright",c("R1","R2","R3"),lty=c(lt1,lt2,lt3),
           col=c(col1,col2,col3), bty="n")


###################################################
### code chunk number 11: GeneralModel.Rnw:219-220
###################################################
Y_r=getAccumulatedRelease(mod) 


###################################################
### code chunk number 12: GeneralModel.Rnw:224-232
###################################################
    plot(t,Y_r[,1],type="l",lty=lt1,col=col1,
    ylab="Accumulated Release (arbitrary
    units)", xlab="Time",
    ylim=c(min(Y_r),max(Y_r)))
    lines(t,Y_r[,2],lt2,type="l",lty=lt2,col=col2)
    lines(t,Y_r[,3],type="l",lty=lt3,col=col3)
    legend("topleft",c("R1","R2","R3"),lty=c(lt1,lt2,lt3),
           col=c(col1,col2,col3), bty="n")


###################################################
### code chunk number 13: GeneralModel.Rnw:269-276
###################################################
Temp=function(t0){ #Temperature in Celsius 
    T0=10   #anual average temperature in Celsius degree
    A=10    #Amplitude in K
    P=1     #Period in years
    T0+A*sin(2*pi*P*t0)
}
plot(t,Temp(t),xlab="Time",ylab="Temperature (degrees Celcius)",type="l")


###################################################
### code chunk number 14: GeneralModel.Rnw:280-288
###################################################
Moist=function(t0){#Moisture in percent
    W0=70       #average moisture in percent
    A=10        #Amplitude of change
    P=1         #Period in years
    ps=pi/7     #phase shift
    W0+A*sin(2*pi*P*t0-ps)
}
plot(t,Moist(t),xlab="Time",ylab="Moisture (percentage)",type="l")


###################################################
### code chunk number 15: GeneralModel.Rnw:301-304
###################################################
xi=function(t0){
    fT.Daycent1(Temp(t0))*as.numeric(fW.Daycent2(Moist(t0))["fRWC"])
}


###################################################
### code chunk number 16: <
###################################################
A_0=matrix(nrow=n,ncol=n,byrow=TRUE,
  c(-0.2,    0,    0, 
     0.1, -0.7,    0,   
     0,    1/2,   -0.5)
)
A_t=BoundLinDecompOp(
  function(t0){xi(t0)*A_0},
  t_start,
  t_end
  )


###################################################
### code chunk number 17: <
###################################################
inputFluxes=function(t0){
    t_peak1=0.75
    t_peak2=1.75
    matrix(nrow=3,
           ncol=1,
           c(
             exp(-((t0-t_peak1)*40.0)^2)+exp(-((t0-t_peak2)*40.0)^2),
             0,
             0
           )
    )
}
inputFluxes_tm=BoundInFlux(
    inputFluxes,
    t_start,
    t_end
) 


###################################################
### code chunk number 18: GeneralModel.Rnw:341-344
###################################################
ifl_1=matrix(nrow=1,ncol=length(t))
for (i in 1:length(t)){ifl_1[i]=inputFluxes(t[i])[1]}
plot(t,ifl_1,xlab="Time",ylab="external inputflux to pool 1",type="l")


###################################################
### code chunk number 19: GeneralModel.Rnw:349-367
###################################################
mod=GeneralModel(t,A_t,c0,inputFluxes_tm)
Y_c=getC(mod)
      plot(t,Y_c[,1],type="l",lty=lt1,col=col1,
           ylab="C stocks (arbitrary units)",
           xlab="Time",
           ylim=c(min(Y_c),1.1*max(Y_c))
      ) 
      lines(t,Y_c[,2],type="l",lty=lt2,col=col2) 
      lines(t,Y_c[,3],type="l",lty=lt3,col=col3) 
      legend(
         "topright",
         c("C in pool 1",
           "C in pool 2",
           "C in pool 3"
         ),
         lty=c(lt1,lt2,lt3),
         col=c(col1,col2,col3), bty="n"
      )


###################################################
### code chunk number 20: GeneralModel.Rnw:393-398
###################################################
    fn="inputFluxForVignetteGeneralModel"
    fn2="inputFluxForVignetteGeneralModelShort"
    subdir=file.path(system.file(package="SoilR"),"extdata")
    p=file.path(subdir,fn)
    p2=file.path(subdir,fn2)


###################################################
### code chunk number 21: GeneralModel.Rnw:401-402
###################################################
#file.show(p)


###################################################
### code chunk number 22: GeneralModel.Rnw:404-430
###################################################
# 1.)
# To make it simpler to keep consistency between the
# data and the vignette 
# we provide the code to create the datafiles here
# To this end we  now create a dataframe using the same time 
# range as all the examples before 
t_peak1=0.75
t_peak2=1.75
df=data.frame("time"=t,"inputFlux"=exp(-((t-t_peak1)*40.0)^2)+exp(-((t-t_peak2)*40.0)^2))

# Additionally we create a second dataset spanning a smaller 
# time range which we will use later to demonstrate 
# the safety net provided by the use of the classes.
d=t_end-t_start
ts2=t_start+d/4
te2=t_end-d/4
t2=seq(ts2,te2,timestep) 
i2=exp(-((t2-t_peak1)*40.0)^2)+exp(-((t2-t_peak2)*40.0)^2)
df2=data.frame(t2,i2)
# temporary uncomment the next lines to write the file to the appropriate source dir and
# comment it out again before you check in the vignette because it will not pass the package check otherwise
#mmdir="~/SoilR/RPackages/SoilR/pkg/inst/extdata"
#write.csv(df ,row.names=FALSE,file.path(mmdir,fn))
#write.csv(df ,row.names=FALSE,p)
#write.csv(df2,row.names=FALSE,file.path(mmdir,fn2))
#write.csv(df2,row.names=FALSE,p2)


###################################################
### code chunk number 23: GeneralModel.Rnw:434-436
###################################################
dfr=read.csv(p)
iTm=BoundInFlux(dfr)


###################################################
### code chunk number 24: GeneralModel.Rnw:440-441
###################################################
mod=GeneralModel(t,A_t,c0,iTm)


###################################################
### code chunk number 25: GeneralModel.Rnw:449-451
###################################################
dfr2=read.csv(p2)
iTm=BoundInFlux(dfr2)


###################################################
### code chunk number 26: GeneralModel.Rnw:455-456
###################################################
#mod=GeneralModel(t,A_t,c0,iTm)


###################################################
### code chunk number 27: GeneralModel.Rnw:464-467
###################################################
getTimeRange(iTm)
min(t)
max(t)


###################################################
### code chunk number 28: GeneralModel.Rnw:478-482
###################################################
ts2
te2
t=seq(ts2,te2,timestep) 
mod=GeneralModel(t,A_t,c0,iTm)


