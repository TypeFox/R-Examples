### R code from vignette source 'PSM.Rnw'

###################################################
### code chunk number 1: Sweavepath
###################################################
sweavepath <- gsub("[\\]",'/',R.home("share/texmf/"))


###################################################
### code chunk number 2: PackageDescription
###################################################
options(keep.source = TRUE, width = 60)
foo <- packageDescription("PSM")
library("PSM")
foo[['Built']] <- gsub("[#$%^&_{}~]",replacement="-", x=foo[['Built']])
options(digits=4)


###################################################
### code chunk number 3: PSMtemplateLinear
###################################################
PSM.template(Linear=TRUE,dimX=3,dimY=1,dimU=0,dimEta=3)


###################################################
### code chunk number 4: DataObj
###################################################
MyData <- list()
MyData[[1]] <- list(Time=1:4,Y=matrix(c(2.1,3.2,3.4,3.7),nrow=1),covar=c(BMI=20.1))
MyData[[2]] <- list(Time=3:7,Y=matrix(c(1.9,2.1,2.0,2.9,3.5),nrow=1),covar=c(BMI=23.4))


###################################################
### code chunk number 5: Redo.flag
###################################################
Redo = FALSE


###################################################
### code chunk number 6: Dose.ModelDefinition
###################################################
Model.SimDose = list()
Model.SimDose$Matrices = function(phi) {
  V1i <- phi$V1i; V2=phi$V2; CL = phi$CL; CLd = phi$CLd;
  matA <- matrix(c(-(CL+CLd)/V1i ,  CLd/V2 ,
                   CLd/V1i , -CLd/V2 ) ,nrow=2,byrow=T)
  matC <- matrix(c(1/V1i,0),nrow=1)
  list(matA=matA,matC=matC)
}
Model.SimDose$X0 = function(Time=Na,phi,U=Na) {
  matrix(0,nrow=2)
}
Model.SimDose$SIG = function(phi) {
  sig1 <- phi[["sig1"]]
  matrix(c( sig1,0,
           -sig1,0), nrow=2, byrow=T)
}
Model.SimDose$S = function(phi) {
  matrix(phi[["S"]])
}
Model.SimDose$h = function(eta,theta,covar) {
  phi <- theta
  phi$V1i <- theta$V1*exp(eta[1])
  phi
}
Model.SimDose$ModelPar = function(THETA){
  V2 <- 10
  CLd <- 0.1
  list(theta=list(V1 = THETA['V1'],V2=V2,CLd=CLd,CL=THETA['CL'], sig1=THETA['sig1'], S=THETA['S']),
       OMEGA=matrix(THETA['OMEGA1']) )
}
SimDose.THETA <-  c(CL=0.05,V1 = 5, sig1 = 10 , S = 20 , OMEGA1 = .2)


###################################################
### code chunk number 7: SamplingScheme
###################################################
N = 5
SimDose.Data <- vector(mode="list",length=N)
for (i in 1:N) {
  SimDose.Data[[i]]$Time <- seq(from=10,by=10,to=400)
  SimDose.Data[[i]]$Dose <-list(
                                Time = c(30,180),
                                State = c(1, 1),
                                Amount = c(1500,1500)
                                )
}


###################################################
### code chunk number 8: Dose.SimulateData
###################################################
if(Redo) {
  SimDose.Data <- PSM.simulate(Model.SimDose, SimDose.Data, SimDose.THETA, deltaTime=.1)
} else
  load("simdose.RData")


###################################################
### code chunk number 9: figure1
###################################################
PSM.plot(SimDose.Data,indiv=1:2,type=c('Y','longX','eta'))
#par(mfcol=c(3,2),mar = c(2, 4, 2, 2)+.1)
#for(id in 1:2) {
#  plot(SimDose.Data[[id]]$Time , SimDose.Data[[id]]$Y,
#         ylab="Observations", main=paste('Individual ',id,', eta= ',
#                                round(SimDose.Data[[id]]$eta,3),sep=""))
#  for(i in 1:2) {
#    plot(SimDose.Data[[id]]$longTime , SimDose.Data[[id]]$longX[i,],type="l",
#         ylab=paste('State',i))
#    rug(SimDose.Data[[id]]$Time)
#  }
#}


###################################################
### code chunk number 10: Dose.EstimateParameters
###################################################
parA <- list(LB=SimDose.THETA*.5, Init=SimDose.THETA , UB=SimDose.THETA*1.5 )
if(Redo) fitA <- PSM.estimate(Model.SimDose,SimDose.Data,parA,CI=T)
fitA[1:5]
SimDose.THETA


###################################################
### code chunk number 11: Dose.EstimateStates
###################################################
if(Redo)
  out <- PSM.smooth(Model.SimDose, SimDose.Data, fitA$THETA, subsample = 20)
# View the data structure
names(out[[1]])


###################################################
### code chunk number 12: figure2
###################################################
par(mfcol=c(3,2),mar = c(2, 4, 2, 2)+.1)
pos <- c('topright','bottomright')
for(id in 1:2) {
  tmp <- out[[id]]$Ys
  plot(SimDose.Data[[id]]$Time , SimDose.Data[[id]]$Y,
         ylab="Observations", main=paste('Individual',id),ylim=range(tmp))
  legend(x=pos[1], legend=c('Smooth est.'),
         lty=c('71A1'),col='blue', xjust=1,yjust=1,bty="n")
  lines(out[[id]]$Time , out[[id]]$Ys,lty='71A1',col='blue')
  for(i in 1:2) {
    plot(SimDose.Data[[id]]$longTime , SimDose.Data[[id]]$longX[i,],type="l",
         ylab=paste('State',i))
    lines(out[[id]]$Time , out[[id]]$Xs[i,],lty='71A1',col="blue")
    rug(SimDose.Data[[id]]$Time)
    legend(x=pos[i], legend=c('Simulation','Smooth est.'),
           lty=c('solid','71A1'),col=c("black","blue"),xjust=1,
           yjust=1,bty="n",y.intersp=.8)
  }
}


###################################################
### code chunk number 13: ISR.ModelDefinition
###################################################
k1 = 0.053; k2 = 0.051; ke = 0.062;
Model.SimISR <- list()
Model.SimISR$Matrices = function(phi) {
  a1  <- phi[["a1"]]
  a2  <- phi[["a2"]]
  B  <- phi[["B"]]
  K  <- phi[["K"]]
  matA <- matrix( c(-(k1+ke) ,  k2 ,   1 ,   0,
                    k1 , -k2 ,   0 ,   0,
                    0 ,   0 , -a1 ,  a1,
                    0 ,   0 ,   0 , -a2),nrow=4,byrow=T)
  matB <- matrix( c(0 , 0   ,
                    0 , 0   ,
                    B , 0   ,
                    0 , a2*K),byrow=T,nrow=4)
  matC <- matrix(c(1,0,0,0),nrow=1)
  matD <- matrix(c(0,0),nrow=1)
  list(matA=matA,matB=matB,matC=matC,matD=matD)
}
Model.SimISR$X0 = function(Time=NA,phi,U=NA) {
  C0 <- phi[["C0"]]
  tmp    <- C0
  tmp[2] <- C0*k1/k2
  tmp[3] <- C0*ke
  tmp[4] <- 0
  matrix(tmp,ncol=1)
}
Model.SimISR$SIG = function(phi) {
  diag( c(0,0,phi[["SIG33"]],0))
}
Model.SimISR$S = function(phi) {
  return( matrix(phi[["S"]]))
}
Model.SimISR$h = function(eta,theta,covar) {
  phi <- theta
  phi[["B"]] <- theta[["B"]]*exp(eta[1])
  phi[["K"]] <- theta[["K"]]*exp(eta[2])
  phi[["C0"]] <- theta[["C0"]]*exp(eta[3])
  return(phi)
}
Model.SimISR$ModelPar = function(THETA){
  list(theta=list(C0=900,S=8500,
                a1=THETA['a1'],a2=THETA['a2'],
                SIG33=THETA['SIG33'],
                K = THETA['K'], B = THETA['B']),
              OMEGA=diag(c(.2,.2,.2))
       )
}



###################################################
### code chunk number 14: ISR.SimulateData
###################################################
Sim.Data <- vector(mode="list",length=2)
for (i in 1:2) {
  Sim.Data[[i]]$Time <- c( 0,15,30,45,60,75,90,120,150,180,210,240,270,300,330,
                          360,420,480,600,615,630,645,660,675,690,720,750,780,810,
                          840,960,1140,1320,1410,1440)
  Sim.Data[[i]]$U <- matrix(c( rep(1,35) ,
                              as.numeric( Sim.Data[[i]]$Time %in% c(0,15,240,600,615)) )
                            ,byrow=T,nrow=2)
}


###################################################
### code chunk number 15: ISR.SimulateData
###################################################
Sim.THETA <-  c(a1=0.02798, a2=0.01048, SIG33=4 , K=427.63 , B=1.7434)
if(Redo) {
  Sim.Data <- PSM.simulate(Model.SimISR, Sim.Data, Sim.THETA, deltaTime=.1 )
} else
  load("simisr.RData")


###################################################
### code chunk number 16: figure3
###################################################
par(mfcol=c(4,2),mar = c(2, 4, 2, 2)+.1)
for(id in 1:2) {
  for(i in 1:4) {
    plot(Sim.Data[[id]]$longTime , Sim.Data[[id]]$longX[i,],type="l", ylab=paste('State',i))
    rug(Sim.Data[[id]]$Time)
    if(i==1) {
      title(main=paste('Individual ',id,', eta=(',paste(round(Sim.Data[[id]]$eta,3),
              collapse=","),')',sep=""))
      points(Sim.Data[[id]]$Time, Sim.Data[[id]]$Y)
      legend(1440, y = max(Sim.Data[[id]]$X[i,]), legend=c('Obs.'),
           pch=21,xjust=1, yjust=1)
    }
  }
}


###################################################
### code chunk number 17: ISR.EstimationModel
###################################################
Model.Est <- list(
            Matrices=function(phi) { list(
              matA=matrix(c(-(k1+ke),  k2, 1,
                             k1     , -k2, 0,
                             0      ,   0, 0  ),ncol=3,byrow=T),
              matB=NA,
              matC=matrix(c(1,0,0),nrow=1),
              matD=NA )  },
            X0 = function(Time=NA,phi=NA,U=NA) {
              C0 <- phi[["C0"]]
              tmp    <- C0
              tmp[2] <- C0*k1/k2
              tmp[3] <- C0*ke
                    return(matrix(tmp,ncol=1) )} ,
            SIG = function(phi) {
                      return( diag( c(1e-3,1e-3,phi[["SIG33"]])) ) } ,
            S = function(phi) {
                      return( matrix(phi[["S"]])) } ,
            h = function(eta,theta,covar) {
              phi <- theta
              phi[["C0"]] <- theta[["C0"]]*exp(eta[1])
              return(phi) } ,
            ModelPar = function(THETA){
              return(list(theta=list(C0=THETA['C0'],S=THETA['S'],SIG33=THETA['SIG33']),
                          OMEGA=matrix(THETA['OMEGA'])))}
            )


###################################################
### code chunk number 18: ISR.RemoveInputFromData
###################################################
Pop.Data <- Sim.Data
for (i in 1:2)
  Pop.Data[[i]]$U <- NULL


###################################################
### code chunk number 19: ISR.ParameterEstimation
###################################################
par1 <- list(LB   = c(C0= 200, S= 50^2, SIG33= 0, OMEGA=.0 ),
             Init = c(C0=1000, S=100^2, SIG33=10, OMEGA=.25),
             UB   = c(C0=3000, S=150^2, SIG33=15, OMEGA=.50))
if(Redo) obj1 <- PSM.estimate(Model.Est, Pop.Data, par1,CI=T,trace=1)
obj1[1:5]


###################################################
### code chunk number 20: ISR.EstimateStates
###################################################
if(Redo)
  Data.Sm <- PSM.smooth( Model.Est , Pop.Data, obj1$THETA, subsample=10)


###################################################
### code chunk number 21: figure4
###################################################
par(mfcol=c(4,2),mar = c(2, 4, 2, 2)+.1)
for(id in 1:2) {
  for(i in 1:4) {
    plot(Sim.Data[[id]]$longTime , Sim.Data[[id]]$longX[i,],type="l",
         ylab=paste('State',i), xlab=paste('Individual',id))
    if(i<=3) lines(Data.Sm[[id]]$Time, Data.Sm[[id]]$Xs[i,],col="blue",lty='71A1')
    rug(Sim.Data[[id]]$Time)
    legend(1440, y = max(Sim.Data[[id]]$X[i,]), legend=c('Simul.','Sm. est.'),
       lty=c('solid','71A1'),col=c("black","blue"),xjust=1,
       yjust=1,bty="n",y.intersp=.8)
    if(i==1)  {
      title(main=paste('Individual',id))
      points(Sim.Data[[id]]$Time, Sim.Data[[id]]$Y)
    }
  }
}


###################################################
### code chunk number 22: figure5
###################################################
par(mfcol=c(2,1),mar = c(3,3, 0, 2)+.1)
for(ID in 1:2) {
  CMT <- 3
  Data <- Data.Sm[[ID]]
  polyCI = c(Data$Xs[CMT,]+sqrt(abs(Data$Ps[CMT,CMT,])),
    rev(Data$Xs[CMT,]-sqrt(abs(Data$Ps[CMT,CMT,]))))
  ymax=max(polyCI)
  plot.new()
  plot.window(xlim=range(Data$Time),ylim=c(0,ymax))
  axis(1);axis(2);box()
  title(xlab="minutes",ylab="pmol/min",line=2)
  polygon( c(Data$Time,rev(Data$Time)) ,polyCI,col=rgb(.9,.9,1),border=NA)
  polygon( c(Data$Time,rev(Data$Time)) ,polyCI,col=rgb(.6,.6,1),density=0)
  lines( Data$Time, Data$Xs[CMT,], type="l",lwd=2,col="blue")
  lines( Sim.Data[[ID]]$longTime, Sim.Data[[ID]]$longX[3,],lwd=.5)
  #points( Sim.Data[[ID]]$Time, Sim.Data[[ID]]$X[3,],pch=20,col=1,cex=.7)
  legend(1440, y = max(polyCI), legend=c('Simulation','Smooth est.'),
         col=c("black","blue"),xjust=1,
         yjust=1,bty="n",y.intersp=.8,lwd=c(1,2),pch=c(NA,NA))
  box()
}


