### LongData contient des trajectoires de plusieurs variables et leur identifiant
### Une trajectoire est une suite de numeric.
### Les trajectoires pour n variables est un array en 3D :
###  - Chaque plan horizontal (premier dimension) est n trajectoires pour un individu
###  - Chaque plan vertical de coupe (deuxième dimension) est un temps Ti
###  - Chaque plan vertical longitudinal (troisième dimension) est une variable
### Les autres variables
###  - idAll est l'ensemble de tous les identifiants
###  - idFewNA est l'ensemble des identifiants des trajectoires non exclues pour cause de manquantes
###  - time est le temps ou les mesures ont été faite.
###  - varName est les noms des variables.
###  - dimTraj est la dimension de la matrice
###  - maxNA est le nombre de NA tolérant avant exclusion de la trajectoire
###  - reverse contiennet les parametres permettant de restaurer les valeurs initiales des trajectoires.


source("testLongData.data.r")


cat("\n####################################################################
########################### Test  LongData #########################
############################# Accesseurs ###########################
####################################################################\n")

cleanProg(.longData.get)
cleanProg(.longData.show)
cleanProg(.longData.print)
cleanProg(.longData.plot)
cleanProg(resizePartition)
cleanProg(calculTrajMean,,,2) # meanNA tapply
cleanProg(calculTrajMeanPoint)
cleanProg(legendCol)
cleanProg(.LongData.Partition.plot,,,1) #LETTERS
cleanProg(varNumAndName)
cleanProg(.LongData.plot3d)
cleanProg(adjustGraph3d)
cleanProg(.LongData.Partition.plot3d,,,1) # LETTERS


ld8["idAll"][1:100]
ld8n["idAll"][1:100]
ld7["idAll"]
ld7n["idAll"]
ld6["idAll"]
ld6n["idAll"]
ld3n["idAll"][1:100]
ld2["idAll"]
ld2n["idAll"]
ld1["idAll"]
ld1n["idAll"]
ld0["idAll"]

ld8["idFewNA"][1:100]
ld7n["idFewNA"]
ld6["idFewNA"]
ld3n["idFewNA"][1:100]
ld2["idFewNA"]
ld1n["idFewNA"]
ld0["idFewNA"]

ld8["varNames"]
ld7n["varNames"]
ld6["varNames"]
ld3n["varNames"]
ld2["varNames"]
ld1n["varNames"]
ld0["varNames"]

ld8["time"]
ld7n["time"]
ld6["time"]
ld3n["time"]
ld2["time"]
ld1n["time"]
ld0["time"]

#ld8["traj"]
#ld7n["traj"]
ld6["traj"]
#ld3n["traj"]
ld2["traj"]
ld1n["traj"]
ld0["traj"]

ld8["dimTraj"]
ld7n["nbVar"]
ld6["nbTime"]
ld3n["nbIdFewNA"]

ld8["maxNA"]
ld7n["maxNA"]
ld6["maxNA"]
ld3n["maxNA"]
ld2["maxNA"]
ld1n["maxNA"]
ld0["maxNA"]

ld8["reverse"]
ld7n["reverse"]
ld6["reverse"]
ld3n["reverse"]
ld2["reverse"]
ld1n["reverse"]
ld0["reverse"]


cat("\n###################################################################
########################## Test  LongData #########################
############################# Affichage ###########################
###################################################################\n")

cleanProg(.longData.show)
ld0
ld1n
ld7
ld8n

cleanProg(.longData.print)
print(ld1)
print(ld2n)
print(ld3n)
print(ld6n)
print(ld7)


cleanProg(resizePartition)
cleanProg(calculTrajMean,,,2)   # meanNA tapply
cleanProg(calculTrajMeanPoint)

cleanProg(varNumAndName)
cleanProg(adjustGraph3d)

cleanProg(.longData.plot)
plot(ld2)
plot(ld3n)
plot(ld1,paramWindows=windowsCut(c(2,1)))
plot(ld1,paramWindows=windowsCut(c(1,2),FALSE))
plot(ld1,paramWindows=windowsCut(c(1,2),TRUE))
plot(ld3n,paramWindows=windowsCut(4))
plot(ld8n,paramTraj=parTraj(col=rep(2:7,20),type="b",pch="4"))
plot(ld8n,paramTraj=parTraj(col=rep(2:7,20),type="b",pch="3"),nbSample=Inf)
plot(ld7n,paramTraj=parTraj(col=rep(2:7,20)))
plot(ld4)
plot(ld8n)


cleanProg(.LongData.plot3d)
plot3d(ld1)
plot3d(ld3n)
plot3d(ld7,paramTraj=parTraj(col=(rep(2:7,20))))
plot3d(ld8n)
#plot3d(ld8n,nbSample=Inf)
plot3d(ld4)

plot3d(ld4,varY=5,varZ=7)
plot3d(ld4,varY="V4")
try(plot3d(ld4,varY=30))
try(plot3d(ld4,varY="VT4"))
plot3d(ld8n,p8cn)

cleanProg(.LongData.Partition.plot,,,1)
plot(ld2,p2a)
plot(ld3n,p3b)
plot(ld1,p1b,paramWindows=windowsCut(c(2,2)))
plot(ld3n,p3a,paramWindows=windowsCut(c(1,3)))
plot(ld8n,p8a,paramTraj=parTraj(col=(rep(2:7,20)),type="b",pch="letters"))
plot(ld8n,p8b,paramTraj=parTraj(col=(rep(2:7,20)),type="b",pch="symbols"),nbSample=Inf)
plot(ld8n,p8b,paramTraj=parTraj(col='clusters',type="b",pch="symbols"),nbSample=Inf)
plot(ld7,p7c,paramTraj=parTraj(col='clusters',type="b",pch="symbols"))
plot(ld7n,p7c,paramTraj=parTraj(col='clusters',type="b",pch="symbols"))
plot(ld7,p7cn,paramTraj=parTraj(col='clusters',type="b",pch="symbols"))
plot(ld7n,p7cn,paramTraj=parTraj(col='clusters',type="b",pch="symbols"))
plot(ld4,p4a)

cleanProg(.LongData.Partition.plot3d,,,1)
plot3d(ld1,p1a)
plot3d(ld3n,p3d)
plot3d(ld3n,p3e)
plot3d(ld7,p7b,paramTraj=parTraj(col=(rep(2:7,20))))
plot3d(ld8n,p8cn)
#plot3d(ld8n,nbSample=Inf)
plot3d(ld4,p4a)

plot3d(ld4,p4b,varY=5,varZ=7)
plot3d(ld4,p4c,varY="V4")
try(plot3d(ld4,p4a,varY=30))
try(plot3d(ld4,p4b,varY="VT4"))
plot3d(ld8n,p8cn)


cat("\n###################################################################
########################## Test  LongData #########################
############################ imputation ###########################
###################################################################\n")

dev.off()
a1 <- c(NA,2,NA,4,NA,1,NA,NA,NA,-5,NA,3,NA,NA)
a2 <- c(NA,NA,4,NA,NA,NA,NA)
a3 <- c(NA,NA,-4,5,NA,NA,NA)
a4 <- c(NA,NA,NA,NA,NA,NA)
a5 <- c(1,NA)
a6 <- c(NA,-1)
a7 <- c(1)

cat("### LOCF: sous fonctions ###\n")

cleanProg(trajImput.LOCB.begin,,,0)
cleanProg(trajImput.LOCF.middle,,,0)

a1a <- trajImput.LOCB.begin(a1)
(a1a <- trajImput.LOCF.middle(a1a))

a2a <- trajImput.LOCB.begin(a2)
(a2a <- trajImput.LOCF.middle(a2a))

a3a <- trajImput.LOCB.begin(a3)
(a3a <- trajImput.LOCF.middle(a3a))

a4a <- trajImput.LOCB.begin(a4)
(a4a <- trajImput.LOCF.middle(a4a))

a5a <- trajImput.LOCB.begin(a5)
(a5a <- trajImput.LOCF.middle(a5a))

a6a <- trajImput.LOCB.begin(a6)
(a6a <- trajImput.LOCF.middle(a6a))

a7a <- trajImput.LOCB.begin(a7)
(a7a <- trajImput.LOCF.middle(a7a))


cat("### LOCF: function complete ###\n")

cleanProg(trajImput.LOCF,,,0)

par(mfrow=c(3,5))
a1aC <- trajImput.LOCF(a1)
plot(a1aC,type="o")
lines(a1,lwd=3,col=2,type="o")

a2aC <- trajImput.LOCF(a2)
plot(a2aC,type="o")
lines(a2,lwd=3,col=2,type="o")

a3aC <- trajImput.LOCF(a3)
plot(a3aC,type="o")
lines(a3,lwd=3,col=2,type="o")

a4aC <- trajImput.LOCF(a4)

a5aC <- trajImput.LOCF(a5)
plot(a5aC,type="o")
lines(a5,lwd=3,col=2,type="o")

a6aC <- trajImput.LOCF(a6)
plot(a6aC,type="o")
lines(a6,lwd=3,col=2,type="o")

trajImput.LOCF(a7)


cat("### LOCB: Sous fonctions ###\n")
cleanProg(trajImput.LOCF.end,,,0)
cleanProg(trajImput.LOCB.middle,,,0)

a1b <- trajImput.LOCF.end(a1)
(a1b <- trajImput.LOCB.middle(a1b))

a2b <- trajImput.LOCF.end(a2)
(a2b <- trajImput.LOCB.middle(a2b))

a3b <- trajImput.LOCF.end(a3)
(a3b <- trajImput.LOCB.middle(a3b))

a4b <- trajImput.LOCF.end(a4)
a4b <- trajImput.LOCB.middle(a4b)

a5b <- trajImput.LOCF.end(a5)
(a5b <- trajImput.LOCB.middle(a5b))

a6b <- trajImput.LOCF.end(a6)
(a6b <- trajImput.LOCB.middle(a6b))

a7b <- trajImput.LOCF.end(a7)
(a7b <- trajImput.LOCB.middle(a7b))


cat("### LOCB: function complete ###\n")
cleanProg(trajImput.LOCB,,,0)

a1bC <- trajImput.LOCB(a1)
plot(a1bC,type="o")
lines(a1,lwd=3,col=2,type="o")

a2bC <- trajImput.LOCB(a2)
plot(a2bC,type="o")
lines(a2,lwd=3,col=2,type="o")

a3bC <- trajImput.LOCB(a3)
plot(a3bC,type="o")
lines(a3,lwd=3,col=2,type="o")

a4bC <- trajImput.LOCB(a4)

a5bC <- trajImput.LOCB(a5)
plot(a5bC,type="o")
lines(a5,lwd=3,col=2,type="o")

a6bC <- trajImput.LOCB(a6)
plot(a6bC,type="o")
lines(a6,lwd=3,col=2,type="o")

trajImput.LOCB(a7)



cat("### linear interpolation 2 : global slope, sous fonctions ###\n")
cleanProg(trajImput.interpoLin.middle,,,0)
cleanProg(trajImput.globalSlope.beginEnd,,,0)

a1c <- trajImput.globalSlope.beginEnd(a1)
(a1c <- trajImput.interpoLin.middle(a1c))

a2c <- trajImput.globalSlope.beginEnd(a2)
a2c <- trajImput.interpoLin.middle(a2c)

a3c <- trajImput.globalSlope.beginEnd(a3)
(a3c <- trajImput.interpoLin.middle(a3c))

a4c <- trajImput.globalSlope.beginEnd(a4)
try(a4c <- trajImput.interpoLin.middle(a4c))

a5c <- trajImput.globalSlope.beginEnd(a5)
a5c <- trajImput.interpoLin.middle(a5c)

a6c <- trajImput.globalSlope.beginEnd(a6)
a6c <- trajImput.interpoLin.middle(a6c)

a7c <- trajImput.globalSlope.beginEnd(a7)
(a7c <- trajImput.interpoLin.middle(a7c))


cat("### linear interpolation 2 : global slope, fonction complete ###\n")
cleanProg(trajImput.linInterGlobal,,,0)

a1cC <- trajImput.linInterGlobal(a1)
plot(a1cC,type="o")
lines(a1,lwd=3,col=2,type="o")

a2cC <- trajImput.linInterGlobal(a2)
plot(a2cC,type="o")
lines(a2,lwd=3,col=2,type="o")

a3cC <- trajImput.linInterGlobal(a3)
plot(a3cC,type="o")
lines(a3,lwd=3,col=2,type="o")

trajImput.linInterGlobal(a4)

a5cC <- trajImput.linInterGlobal(a5)
plot(a5cC,type="o")
lines(a5,lwd=3,col=2,type="o")

a6cC <- trajImput.linInterGlobal(a6)
plot(a6cC,type="o")
lines(a6,lwd=3,col=2,type="o")

a7cC <- trajImput.linInterGlobal(a7)



cat("### Linear interpolation 3 : Local slope, sous fonctions ###\n")
cleanProg(trajImput.localSlope.beginEnd,,,0)

a1d <- trajImput.localSlope.beginEnd(a1)
(a1d <- trajImput.interpoLin.middle(a1d))

a2d <- trajImput.localSlope.beginEnd(a2)
try(a2d <- trajImput.interpoLin.middle(a2d))

a3d <- trajImput.localSlope.beginEnd(a3)
(a3d <- trajImput.interpoLin.middle(a3d))

a4d <- trajImput.localSlope.beginEnd(a4)
try(a4d <- trajImput.interpoLin.middle(a4d))

a5d <- trajImput.localSlope.beginEnd(a5)
try(a5d <- trajImput.interpoLin.middle(a5d))

a6d <- trajImput.localSlope.beginEnd(a6)
try(a6d <- trajImput.interpoLin.middle(a6d))

a7d <- trajImput.localSlope.beginEnd(a7)
(a7d <- trajImput.interpoLin.middle(a7d))


cat("###  Linear interpolation 3 : Local slope, fonction complete ###\n")
cleanProg(trajImput.linInterLocal,,,0)

a1dC <- trajImput.linInterLocal(a1)
plot(a1dC,type="o")
lines(a1,lwd=3,col=2,type="o")

a2dC <- trajImput.linInterLocal(a2)
plot(a2dC,type="o")
lines(a2,lwd=3,col=2,type="o")

a3dC <- trajImput.linInterLocal(a3)
plot(a3dC,type="o")
lines(a3,lwd=3,col=2,type="o")

try(trajImput.linInterLocal(a4))

a5dC <- trajImput.linInterLocal(a5)
plot(a5dC,type="o")
lines(a5,lwd=3,col=2,type="o")

a6dC <- trajImput.linInterLocal(a6)
plot(a6dC,type="o")
lines(a6,lwd=3,col=2,type="o")

a7dC <- trajImput.linInterLocal(a7)



cat("### Linear interpolation 4 : LOCF, fonction complete ###\n")
cleanProg(trajImput.linInterLOCF,,,0)

a1eC <- trajImput.linInterLOCF(a1)
plot(a1eC,type="o")
lines(a1,lwd=3,col=2,type="o")

a2eC <- trajImput.linInterLOCF(a2)
plot(a2eC,type="o")
lines(a2,lwd=3,col=2,type="o")

a3eC <- trajImput.linInterLOCF(a3)
plot(a3eC,type="o")
lines(a3,lwd=3,col=2,type="o")

try(trajImput.interpoLin(a4))

a5eC <- trajImput.linInterLOCF(a5)
plot(a5eC,type="o")
lines(a5,lwd=3,col=2,type="o")

a6eC <- trajImput.linInterLOCF(a6)
plot(a6eC,type="o")
lines(a6,lwd=3,col=2,type="o")

a7eC <- trajImput.linInterLOCF(a7)



cat("### Linear interpolation : bissectrice, sous fonctions ###\n")
cleanProg(trajImput.bissectrice.beginEnd,,,0)

a1f <- trajImput.bissectrice.beginEnd(a1)
(a1f <- trajImput.interpoLin.middle(a1f))

a2f <- trajImput.bissectrice.beginEnd(a2)
(a2f <- trajImput.interpoLin.middle(a2f))

a3f <- trajImput.bissectrice.beginEnd(a3)
(a3f <- trajImput.interpoLin.middle(a3f))

a4f <- trajImput.bissectrice.beginEnd(a4)
try(a4f <- trajImput.interpoLin.middle(a4f))

a5f <- trajImput.bissectrice.beginEnd(a5)
a5f <- trajImput.interpoLin.middle(a5f)

a6f <- trajImput.bissectrice.beginEnd(a6)
a6f <- trajImput.interpoLin.middle(a6f)

a7f <- trajImput.bissectrice.beginEnd(a7)
(a7f <- trajImput.interpoLin.middle(a7f))


cat("###  Linear interpolation : bissectrice, fonction complete ###\n")
cleanProg(trajImput.linInterBissectrice,,,0)
cleanProg(bissectrice)

a1fC <- trajImput.linInterBissectrice(a1)
plot(a1fC,type="o")
lines(a1,lwd=3,col=2,type="o")

a2fC <- trajImput.linInterBissectrice(a2)
plot(a2fC,type="o")
lines(a2,lwd=3,col=2,type="o")

a3fC <- trajImput.linInterBissectrice(a3)
plot(a3fC,type="o")
lines(a3,lwd=3,col=2,type="o")

try(trajImput.linInterBissectrice(a4))

a5fC <- trajImput.linInterBissectrice(a5)
plot(a5fC,type="o")
lines(a5,lwd=3,col=2,type="o")

a6fC <- trajImput.linInterBissectrice(a6)
plot(a6fC,type="o")
lines(a6,lwd=3,col=2,type="o")

a7fC <- trajImput.linInterBissectrice(a7)


cat("###  Imputation, fonction complete sur array ###\n")
cleanProg(trajImputArray,,,7) # rangeNA trajImput.linInterBissectrice trajImput.linInterGlobal trajImput.linInterLocal trajImput.linInterLOCF trajImput.LOCB trajImput.LOCF
cleanProg(.trajImputLongData,,,1) #rangeNA
(tr1n)
imputation(tr1n,method="LOCF")
imputation(tr1n,method="LOCB")
imputation(tr1n,method="LI-Bissectrice")
(tr1n)
imputation(tr1n,method="LI-Global")
imputation(tr1n,method="LI-Local")
imputation(tr1n,method="LI-LOCBF")

(tr3n)
imputation(tr3n,method="LOCF")
imputation(tr3n,method="LOCB")
imputation(tr3n,method="LI-Bissectrice")
(tr3n)
imputation(tr3n,method="LI-Global")
imputation(tr3n,method="LI-Local")
imputation(tr3n,method="LI-LOCBF")

(tr4n)
imputation(tr4n,method="LOCF")
imputation(tr4n,method="LOCB")
imputation(tr4n,method="LI-Bissectrice")
(tr4n)
imputation(tr4n,method="LI-Global")
imputation(tr4n,method="LI-Local")
imputation(tr4n,method="LI-LOCBF")


imputation(tra1,method="LOCF")
imputation(tra2,method="LOCB")
imputation(tra3,method="LI-Bissectrice")
imputation(tra1,method="LI-Global")
imputation(tra2,method="LI-Local")
imputation(tra3,method="LI-LOCBF")


cat("###  Imputation, fonction complete sur LongData ###\n")

(ldA<-ldB<-ldC<-ld1n)
imputation(ldA,method="LOCF")
imputation(ldB,method="LOCB")
imputation(ldC,method="LI-Bissectrice")
(ldA<-ldB<-ldC<-ld1n)
imputation(ldA,method="LI-Global")
imputation(ldB,method="LI-Local")
imputation(ldC,method="LI-LOCBF")

(ldA<-ldB<-ldC<-ld2n)
imputation(ldA,method="LOCF")
imputation(ldB,method="LOCB")
imputation(ldC,method="LI-Bissectrice")
(ldA<-ldB<-ldC<-ld2n)
imputation(ldA,method="LI-Global")
imputation(ldB,method="LI-Local")
imputation(ldC,method="LI-LOCBF")
(ldA<-ldB<-ldC<-ld3n)
imputation(ldA,method="LOCF")
imputation(ldB,method="LOCB")
imputation(ldC,method="LI-Bissectrice")
(ldA<-ldB<-ldC<-ld3n)
imputation(ldA,method="LI-Global")
imputation(ldB,method="LI-Local")
imputation(ldC,method="LI-LOCBF")

imputation(tra1,method="LOCF")
imputation(tra2,method="LOCB")
imputation(tra3,method="LI-Bissectrice")
imputation(tra1,method="LI-Global")
imputation(tra2,method="LI-Local")
imputation(tra3,method="LI-LOCBF")
imputation(tra2,method="LOCF")
imputation(tra3,method="LOCB")
imputation(tra1,method="LI-Bissectrice")
imputation(tra2,method="LI-Global")
imputation(tra3,method="LI-Local")
imputation(tra1,method="LI-LOCBF")


cat("###  Scale ###\n")
cleanProg(.longData.scale,,,2) # meanNA sdNA
scale(ld1)
scale(ld1n)
scale(ld2)
scale(ld2n)
scale(ld3n)
scale(ld4)
scale(ld6)
scale(ld6n)
scale(ld7)
scale(ld7n)
scale(ld8)
scale(ld8n)

ld3n
ld1
ld2

cat("###  restaureRealData ###\n")
cleanProg(.longData.restaureRealData)
restaureRealData(ld1)
restaureRealData(ld1n)
restaureRealData(ld2)
restaureRealData(ld2n)
restaureRealData(ld3n)
restaureRealData(ld4)
restaureRealData(ld6)
restaureRealData(ld6n)
restaureRealData(ld7)
restaureRealData(ld7n)
restaureRealData(ld8)
restaureRealData(ld8n)

ld1
ld2



cat("\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
++++++++++++++++++++++++ Fin Test  LongData ++++++++++++++++++++++++
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")
