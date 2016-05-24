### R code from vignette source 'NAPPA.Snw'

###################################################
### code chunk number 1: NAPPA.Snw:21-26
###################################################

require(NAPPA)
data(NS.Lung)

NAPPA.Lung.simple <- NAPPA(NS.Lung , tissueType = "tumour")


###################################################
### code chunk number 2: NAPPA.Snw:31-32
###################################################
lungsamplecolours <- sample(colors(),ncol(NAPPA.Lung.simple))


###################################################
### code chunk number 3: fig1
###################################################
plot(prcomp(t(NAPPA.Lung.simple))$x[,1:2], pch=21 , 
bg=lungsamplecolours , cex=1.5)


###################################################
### code chunk number 4: NAPPA.Snw:45-47
###################################################
NAPPA.Lung.detailed <- NAPPA(NS.Lung , tissueType = "tumour" ,
 output=c("All","Steps"))


###################################################
### code chunk number 5: NAPPA.Snw:52-53
###################################################
identical( NAPPA.Lung.simple , NAPPA.Lung.detailed$GeneExpression )


###################################################
### code chunk number 6: fig2
###################################################
plot(NAPPA.Lung.detailed$PosFactor , NAPPA.Lung.detailed$HousekeepingFactor, 
xlab="Positive Control Factor" , ylab="Housekeeping Normalisation Factor" , 
pch=21 , bg=lungsamplecolours  , cex=1.5)


###################################################
### code chunk number 7: NAPPA.Snw:72-74
###################################################
NAPPA.Lung.batch1 <- NAPPA(NS.Lung[,1:5] , tissueType = "tumour" , 
output="All")


###################################################
### code chunk number 8: fig3
###################################################
plot(NAPPA.Lung.batch1$GeneExpression , NAPPA.Lung.simple[,1:2],
xlab="Interim Results" , ylab="Final Results" , main="First 2 Samples")
abline(0 , 1 , col="RED" , lwd=2)


###################################################
### code chunk number 9: NAPPA.Snw:87-88
###################################################
identical( NAPPA.Lung.batch1$GeneExpression , NAPPA.Lung.simple[,1:2] )


###################################################
### code chunk number 10: NAPPA.Snw:93-94
###################################################
NAPPA.Lung.frozen1 <- NAPPA(NS.Lung , tissueType = "tumour" , NReferenceSamples = 2)


###################################################
### code chunk number 11: fig4
###################################################
plot(NAPPA.Lung.batch1$GeneExpression , NAPPA.Lung.frozen1[,1:2], 
xlab="Interim Results", ylab="Final Results Based On Frozen Analysis", 
main="First 2 Samples")
abline(0 , 1 , col="RED" , lwd=2)


###################################################
### code chunk number 12: NAPPA.Snw:106-107
###################################################
identical( NAPPA.Lung.batch1$GeneExpression , NAPPA.Lung.frozen1[,1:2] )


###################################################
### code chunk number 13: NAPPA.Snw:112-115
###################################################
NAPPA.Lung.frozen2 <- NAPPA(NS.Lung , tissueType = "tumour" ,
 betas=NAPPA.Lung.batch1$Betas ,
 hknormfactor.mean=NAPPA.Lung.batch1$HousekeepingFactor.Mean)


###################################################
### code chunk number 14: fig5
###################################################
plot(NAPPA.Lung.batch1$GeneExpression , NAPPA.Lung.frozen2[,1:2], 
xlab="Interim Results" , ylab="Final Results Based On Frozen Analysis", 
main="First 2 Samples")
abline(0 , 1 , col="RED" , lwd=2)


###################################################
### code chunk number 15: NAPPA.Snw:127-128
###################################################
identical( NAPPA.Lung.batch1$GeneExpression , NAPPA.Lung.frozen2[,1:2] )


###################################################
### code chunk number 16: NAPPA.Snw:136-139
###################################################
data(NS.Dilution)
NAPPA.Dilution <- NAPPA(NS.Dilution , betas=NAPPA.Lung.detailed$Betas, 
output="All")


###################################################
### code chunk number 17: NAPPA.Snw:144-146
###################################################
dilutioncolours <- rep(heat.colors(6),4)
dilutionshapes <- rep(c(21,22,24,25) , each=6)


###################################################
### code chunk number 18: fig6
###################################################
plot(prcomp(t(NAPPA.Dilution$GeneExpression))$x[,1:2], 
pch=dilutionshapes , bg=dilutioncolours , cex=1.5)


###################################################
### code chunk number 19: fig7
###################################################
opar <- par(mfrow=c(2,1) , mar=c(1,3,4,1))
plot(c(1,6) , range(NAPPA.Dilution$HousekeepingFactor) , type="n", 
xlab="" , ylab="" , axes=F , main="Housekeeping Normalisation Factors")
box()
axis(2 , las=2)
for(i in 1:4) {
  lines(1:6 , NAPPA.Dilution$HousekeepingFactor[6*i + (-5:0)] , col=i , lwd=2)
  points(1:6 , NAPPA.Dilution$HousekeepingFactor[6*i + (-5:0)] , bg=i , pch=21)
}


###################################################
### code chunk number 20: fig8
###################################################
plot(c(1,6) , range(NAPPA.Dilution$PosFactor) , type="n", 
xlab="" , ylab="" , axes=F , main="Positive Control Factors")
box()
axis(2 , las=2)
for(i in 1:4) {
  lines(1:6 , NAPPA.Dilution$PosFactor[6*i + (-5:0)] , col=i , lwd=2)
  points(1:6 , NAPPA.Dilution$PosFactor[6*i + (-5:0)] , bg=i , pch=21)
}
par(opar)


###################################################
### code chunk number 21: NAPPA.Snw:191-193
###################################################
NAPPA.Dilution.subtract <- NAPPA(NS.Dilution , hk.method="subtract", 
output="All")


###################################################
### code chunk number 22: NAPPA.Snw:198-205
###################################################
blandaltmanaplot <- function(x,y,...) {
  d <- x-y
  a <- 0.5*(x+y)
  plot( a , d , pch="." , ...)
  abline(h=0 , col="BLUE")
  lines(smooth.spline(y=d , x=a , df=4) , col="RED")
}


###################################################
### code chunk number 23: fig9
###################################################
opar <- par(mfcol=c(2,5) , mar=c(4,4,2,0) , oma=c(0,0,0.5,0.5))
for(i in 1:5) 
{
blandaltmanaplot(c(NAPPA.Dilution.subtract$GeneExpression[,i+c(1,7,13,19)]),
c(NAPPA.Dilution.subtract$GeneExpression[,c(1,7,13,19)]) , ylim=c(-5,10), 
xlim=c(0,15) , ylab="" , xlab="") 
title(paste(c(50,25,12.5,6.25,3.125)[i],"ng",collapse="") , line=1, 
cex.main=1.5)
if(i==1) title(ylab="Difference From 100ng Expression\nSubtraction Normalisation" ,
line=2)

blandaltmanaplot( c(NAPPA.Dilution$GeneExpression[,i+c(1,7,13,19)]), 
c(NAPPA.Dilution$GeneExpression[,c(1,7,13,19)]) , ylim=c(-5,10),
xlim=c(0,15) , ylab="" , xlab="") 
if(i==1) title(ylab="Difference From 100ng Expression\nSubtraction Normalisation" ,
line=2)
if(i==3) title(xlab="Average Of 100ng &\nDilution Expression")
}
par(opar)


###################################################
### code chunk number 24: NAPPA.Snw:234-235
###################################################
data(nsolver.Dilution)


###################################################
### code chunk number 25: NAPPA.Snw:240-254
###################################################
require(NanoStringNorm)
NS.Dilution.NSN <- cbind(NS.Dilution[-(1:21),1:3], 
sapply(NS.Dilution[-(1:21),-(1:3)] , as.numeric))
colnames(NS.Dilution.NSN)[1:3] <- c('Code.Class', 'Name', 'Accession')
colnames(NS.Dilution.NSN)[4:27] <- colnames(nsolver.Dilution)
NSN.Dilution <- NanoStringNorm(NS.Dilution.NSN , CodeCount = 'geo.mean',
                           Background = 'mean',
                           SampleContent = 'housekeeping.geo.mean',
                           round.values = TRUE,
                           take.log = TRUE)$normalized.data
NSN.Dilution <- NSN.Dilution[NSN.Dilution$Code.Class=="Endogenous",-(1:3)]

all(rownames(NSN.Dilution) == rownames(NAPPA.Dilution$GeneExpression))
nsolver.Dilution <- nsolver.Dilution[rownames(NAPPA.Dilution$GeneExpression),]


###################################################
### code chunk number 26: NAPPA.Snw:259-267
###################################################
NAPPA.Dilution.NSN <-  NAPPA(NS.Dilution , scaleFOV=F,
background.method="subtract.global" , nposcontrols=6 , 
poscontrol.method="geometric.mean" , hk.method="subtract",
output="All")
NAPPA.Dilution.NS <- NAPPA(NS.Dilution , scaleFOV=F, 
background.method="none" , nposcontrols=6,
poscontrol.method="geometric.mean" , hk.method="subtract",
output="All")


###################################################
### code chunk number 27: fig10
###################################################
opar <- par(mfrow=c(3,2) , mar=c(0,0,2,0))
plot(prcomp(t(NAPPA.Dilution$GeneExpression))$x[,1:2] , pch=dilutionshapes,
bg=dilutioncolours , cex=1.5 , axes=F , main="NAPPA")
box()
plot(prcomp(t(NAPPA.Dilution$GeneExpression))$x[,1:2] , pch=dilutionshapes,
bg=dilutioncolours , cex=1.5 , axes=F , main="NAPPA")
box()
plot(prcomp(t(NSN.Dilution))$x[,1:2] , pch=dilutionshapes , bg=dilutioncolours,
cex=1.5 , axes=F , main="NanoStringNorm")
box()
plot(prcomp(t(NAPPA.Dilution.NSN$GeneExpression))$x[,1:2] , pch=dilutionshapes,
bg=dilutioncolours , cex=1.5 , axes=F , main="NAPPA version of NanoStringNorm")
box()
plot(prcomp(t(log2(nsolver.Dilution)))$x[,1:2] , pch=dilutionshapes,
bg=dilutioncolours , cex=1.5 , axes=F , main="nsolver")
box()
plot(prcomp(t(NAPPA.Dilution.NS$GeneExpression))$x[,1:2] , pch=dilutionshapes,
bg=dilutioncolours , cex=1.5 , axes=F , main="NAPPA version of nsolver")
box()
par(opar)


