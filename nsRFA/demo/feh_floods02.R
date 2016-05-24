# ------------------------------------------------------------------------------------- #
# ------------------------------------------------------------------------------------- #
#                     ANALYSIS OF FEH DATA: CLASSIFICATION-VARIABLES                    #
# ------------------------------------------------------------------------------------- #
# ------------------------------------------------------------------------------------- #


# Data loading:

data(FEH1000)



# Criteria used in the FEH to choose stations for pooling groups:
#  n>7
#  area, saar and bfihost are known
#  urbext<0.025
#  area>0.5 \n

n <- tapply(am[,4],am[,1],length)
urbext <- cd[,"urbext1990"]
area <- cd[,"dtm_area"]
cd696 <- cd[(!is.nan(cd[,"dtm_area"]))&(!is.nan(cd[,"saar"]))&(!is.nan(cd[,"bfihost"]))&(n>7)&(urbext<0.025)&(area>0.5),]

fac <- factor(am[,"number"],levels=cd696[,"number"])
am696 <- am[!is.na(fac),]
nlevels(as.factor(am696[,"number"]))


# FEH classification variables:
#  lnAREA, lnSAAR, BFIHOST.

lcv <- tapply(am696[,4],am696[,1],LCV)
lca <- tapply(am696[,4],am696[,1],LCA)

lnAREA <- log(cd696[,"dtm_area"])
lnSAAR <- log(cd696[,"saar"])
BFIHOST <- cd696[,"bfihost"]

lcvVSclassvarFEH <- lm(lcv ~ lnAREA + lnSAAR + BFIHOST)
lcaVSclassvarFEH <- lm(lca ~ lnAREA + lnSAAR + BFIHOST)

graphics.lm <- function (regr) {
 par(mfrow=c(1,2), cex=1)
  plot(regr$fitted.values,regr$model[,1],pch=".",cex=3,xlab="Originals",ylab="Fitted")
  abline(0,1,lty=3)
  intervals <- predinterval.lm(regr)
  intervals <- intervals[order(intervals[,1]),]
  lines(intervals[,c(1,2)],lty=2)
  lines(intervals[,c(1,3)],lty=2)
  Rsq <- signif(R2.lm(regr),3)
  rmse <- signif(RMSE.lm(regr),3)
  rmsep <- signif(RMSEP(regr$model[,1],regr$fitted.values)*100,3)
  mtext(paste("R2 = ",Rsq[1]),3,-1.5,adj=0.05,cex=1)
  mtext(paste("R2adj = ",Rsq[2]),3,-2.5,adj=0.05,cex=1)
  mtext(paste("RMSE = ",rmse),1,-2.5,adj=0.95,cex=1)
  mtext(paste("RMSEP = ",rmsep,"%"),1,-1.5,adj=0.95,cex=1)

  normplot(regr$residuals,pch=".",cex=3,xlab="Residuals")
 par(mfrow=c(1,1),cex=1)
 title(main=paste(names(regr$model[1]),"~",paste(names(regr$model[-1]), collapse=", ")),cex.main=1,font.main=1)
}

graphics.lm(lcvVSclassvarFEH)
par(ask = interactive())

graphics.lm(lcaVSclassvarFEH)



# Alternative classification variables:

NGRX <- cd696[,"ihdtm_ngr_x"]
NGRY <- cd696[,"ihdtm_ngr_y"]
AREA <- cd696[,"dtm_area"]
SAAR <- cd696[,"saar"]
SPRHOST <- cd696[,"sprhost"]
FARL <- cd696[,"farl"]
RMED1D <- cd696[,"rmed_1d"]
RMED2D <- cd696[,"rmed_2d"]
RMED1H <- cd696[,"rmed_1h"]
SMDBAR <- cd696[,"smdbar"]
PROPWET <- cd696[,"propwet"]
LDP <- cd696[,"ldp"]
DPLBAR <- cd696[,"dplbar"]
ALTBAR <- cd696[,"altbar"]
DPSBAR <- cd696[,"dpsbar"]
ASPBAR <- cd696[,"aspbar"]
ASPVAR <- cd696[,"aspvar"]
URBEXT <- cd696[,"urbext1990"]

potCVnat <- data.frame(cbind(NGRX,NGRY,AREA,SAAR,BFIHOST,SPRHOST,FARL,RMED1D,RMED2D,RMED1H,
                             SMDBAR,PROPWET,LDP,DPLBAR,ALTBAR,DPSBAR,ASPBAR,ASPVAR,URBEXT))
potCVln <- log(potCVnat[-c(17,19)]); names(potCVln) <- paste("ln",names(potCVln),sep="")
potCV <- cbind(potCVnat,potCVln)

# bestregressions <- function(dip,ind) {
#  Y <- as.numeric(dip)
#  X <- ind
#  Sy <- var(Y)
#  Sx <- var(X)
#  Sxy <- var(X,Y)
#  Dm.mat <- Sx
#  Dm.H <- Sxy %*% t(Sxy)/Sy
#  require(subselect)
#  Dm.leaps <- leaps(Dm.mat, kmin=1, kmax=6, H=Dm.H, r=1, nsol=3)
#  Dm.leaps
#  for(i in 6:1) {for(j in 1:3) {print(colnames(X)[Dm.leaps$subsets[j,c(1:6),i]])}}
# }

# bestregressions(lcv,potCVnat[-19])
# bestregressions(lcv,potCVln)
# bestregressions(lcv,cbind(potCVnat[-19],potCVln[-c(1,4:5,7:13)]))
# bestregressions(lcv,cbind(potCVnat[-c(1,4:5,7:13,19)],potCVln))
# bestregressions(lcv,cbind(potCVnat[-c(1:2,19)],potCVln))
names(cbind(potCVnat[-c(1:2,19)],potCVln))

graphics.lm(lm(lcv ~ lnAREA + lnSAAR + lnRMED2D,data=potCV))
graphics.lm(lm(lcv ~ DPLBAR + lnSAAR + lnRMED2D,data=potCV))
graphics.lm(lm(lcv ~ DPLBAR + ALTBAR + lnSAAR + lnRMED2D,data=potCV))
graphics.lm(lm(lcv ~ RMED1D + ALTBAR + lnAREA + lnSAAR + lnRMED2D + lnALTBAR,data=potCV))

print(prt.lm(lm(lcv ~ DPLBAR + lnSAAR + lnRMED2D,data=potCV)))
print(prt.lm(lm(lcv ~ lnAREA + lnSAAR + lnRMED1D + lnRMED2D,data=potCV)))
print(prt.lm(lm(lcv ~ DPLBAR + ALTBAR + lnSAAR + lnRMED2D + lnALTBAR,data=potCV)))
print(prt.lm(lm(lcv ~ RMED1D + ALTBAR + lnAREA + lnSAAR + lnRMED2D + lnALTBAR,data=potCV)))


# Distance matrices and Mantel test approach: best models

# distancesAD <- AD.dist(am696[,4],am696[,1])                                   # 1 hour
# distancesCVnat <- data.frame(apply(potCVnat,2,dist))                          
# distancesCVln <- data.frame(apply(potCVln,2,dist))                            
# distancesCV <- data.frame(apply(cbind(potCVnat[-c(1:2,19)],potCVln),2,dist))  

# bestregressions(as.numeric(distancesAD),distancesCVnat)  
# bestregressions(as.numeric(distancesAD),distancesCVln)   
# bestregressions(as.numeric(distancesAD),distancesCV)     


# graphics.lm(lm(distancesAD ~ FARL + RMED2D + lnSAAR,data=distancesCV))
# graphics.lm(lm(distancesAD ~ FARL + SMDBAR + DPLBAR + lnNGRY + lnSAAR + lnRMED2D,data=distancesCV))
# graphics.lm(lm(distancesAD ~ DPLBAR + lnSAAR + lnRMED2D,data=distancesCV))


# Distance matrices and Mantel test approach: mantel test

# regr1CV <- lm(distancesAD ~ SMDBAR,data=distancesCV)
# regr2CV <- lm(distancesAD ~ RMED2D + lnSAAR,data=distancesCV)
# regr3CV <- lm(distancesAD ~ FARL + RMED2D + lnSAAR,data=distancesCV)
# regr4CV <- lm(distancesAD ~ FARL + RMED2D + DPLBAR + lnSAAR,data=distancesCV)
# regr5CV <- lm(distancesAD ~ FARL + RMED2D + lnNGRY + lnSAAR + lnDPLBAR,data=distancesCV)
# regr6CV <- lm(distancesAD ~ FARL + SMDBAR + DPLBAR + lnNGRY + lnSAAR + lnRMED2D,data=distancesCV)

# mantel.lm(regr1CV, Nperm=100)   # 2 min 
# mantel.lm(regr2CV, Nperm=100)   # 4 min
# mantel.lm(regr3CV, Nperm=100)   # 6 min
# mantel.lm(regr4CV, Nperm=100)   # 6 min
# mantel.lm(regr5CV, Nperm=100)   # 6 min
# mantel.lm(regr6CV, Nperm=100)   # 6 min
# mantel.lm(lm(distancesAD ~ RMED2D + RMED1D,data=distancesCV), Nperm=100)
# prt.lm(lm(lcv ~ FARL + SMDBAR + DPLBAR + lnNGRY + lnSAAR + lnRMED2D,data=potCV))



# Distance matrices and Mantel test approach: regression statistics

# R2.lm(regr1CV)
# R2.lm(regr2CV)
# R2.lm(regr3CV)
# R2.lm(regr4CV)
# R2.lm(regr5CV)
# R2.lm(regr6CV)

# RMSEP(regr1CV$model[,1],regr1CV$fitted.values)*100
# RMSEP(regr2CV$model[,1],regr2CV$fitted.values)*100
# RMSEP(regr3CV$model[,1],regr3CV$fitted.values)*100
# RMSEP(regr4CV$model[,1],regr4CV$fitted.values)*100
# RMSEP(regr5CV$model[,1],regr5CV$fitted.values)*100
# RMSEP(regr6CV$model[,1],regr6CV$fitted.values)*100



# Formation of groups and homogeneity.
# Classification variables:
#  0) lnAREA, lnSAAR, BFIHOST (same weights + lnAREA/sqrt(2))\n

# term1 <- log(cd696[,"dtm_area"])/(sd(log(cd696[,"dtm_area"]))*sqrt(2))
# term2 <- log(cd696[,"saar"])/sd(log(cd696[,"saar"]))
# term3 <- cd696[,"bfihost"]/sd(cd696[,"bfihost"])

# roi.cd <- data.frame(cbind(term1,term2,term3))
# row.names(roi.cd) <- cd696[,"number"]

# roi00.50year <- new.env()
# for(i in 1:696) {
# print(paste(i,"/ 696"))
# assign(as.character(row.names(roi.cd)[i]), roi.st.year(roi.cd[i,],as.data.frame(roi.cd),
#         row.names(roi.cd),am696[,"am"],am696[,"number"],test="HW and AD",station.year=250,Nsim=100), env=roi00.50year)
# }
# roi00.50year <- as.list(roi00.50year)

# estrai.region <- function (x) {x$region}
# estrai.test <- function (x) {x$test}
# regioni.50year.0 <- sapply(roi00.50year, estrai.region)
# test.50year.0 <- sapply(roi00.50year, estrai.test)
# mL.50year.0 <- mean(sapply(regioni.50year.0,length))
# mH1.50year.0 <- mean(test.50year.0["H1",])
# mH2.50year.0 <- mean(test.50year.0["H2",])
# mAD.50year.0 <- median(test.50year.0["P",])
# gH1gr2.50year.0 <- sum(test.50year.0["H1",]>2)/696
# gH1gr4.50year.0 <- sum(test.50year.0["H1",]>4)/696
# gH2gr2.50year.0 <- sum(test.50year.0["H2",]>2)/696
# gH2gr4.50year.0 <- sum(test.50year.0["H2",]>4)/696
# gADgr99.50year.0 <- sum(test.50year.0["P",]>.99)/696
# gADgr95.50year.0 <- sum(test.50year.0["P",]>.95)/696

# table.0 <- signif(c(mL.50year.0,mH1.50year.0,mH2.50year.0,mAD.50year.0,
#            gH1gr2.50year.0*100,gH1gr4.50year.0*100,gH2gr2.50year.0*100,gH2gr4.50year.0*100,
#            gADgr95.50year.0*100,gADgr99.50year.0*100),3)
# names(table.0) <- c("Avg. n sites","m(H1)","m(H2)","med(p(AD))","% H1>2","% H1>4",
#                     "% H2>2","% H2>4","% p(AD)>0.95","% p(AD)>0.99")
# print(table.0)

table.0 <- c(11.20,2.87,1.53,0.99,62.50,25.10,33.80,5.89,68.20,46.00)
names(table.0) <- c("Avg. n sites","m(H1)","m(H2)","med(p(AD))","% H1>2","% H1>4",
                   "% H2>2","% H2>4","% p(AD)>0.95","% p(AD)>0.99")
print(table.0)



# Formation of groups and homogeneity. 
# Classification variables:
#  1) DPLBAR, lnSAAR, lnRMED2D (with the same weight)

# term1 <- cd696[,"dplbar"]/sd(cd696[,"dplbar"])
# term2 <- log(cd696[,"saar"])/sd(log(cd696[,"saar"]))
# term3 <- log(cd696[,"rmed_2d"])/sd(log(cd696[,"rmed_2d"]))

# roi.cd <- data.frame(cbind(term1,term2,term3))
# row.names(roi.cd) <- cd696[,"number"]

# roi.st.year(roi.cd["54088",],roi.cd,row.names(roi.cd),am696[,"am"],
#             am696[,"number"],test="HW and AD",station.year=250,Nsim=500)

# roi01.50year <- new.env()
# for(i in 1:696) {
# print(paste(i,"/ 696"))
# assign(as.character(row.names(roi.cd)[i]), roi.st.year(roi.cd[i,],as.data.frame(roi.cd),
#         row.names(roi.cd),am696[,"am"],am696[,"number"],test="HW and AD",station.year=250,Nsim=100), env=roi01.50year)
# }
# roi01.50year <- as.list(roi01.50year)

# estrai.region <- function (x) {x$region}
# estrai.test <- function (x) {x$test}
# regioni.50year.1 <- sapply(roi01.50year, estrai.region)
# test.50year.1 <- sapply(roi01.50year, estrai.test)
# mL.50year.1 <- mean(sapply(regioni.50year.1,length))
# mH1.50year.1 <- mean(test.50year.1["H1",])
# mH2.50year.1 <- mean(test.50year.1["H2",])
# mAD.50year.1 <- median(test.50year.1["P",])
# gH1gr2.50year.1 <- sum(test.50year.1["H1",]>2)/696
# gH1gr4.50year.1 <- sum(test.50year.1["H1",]>4)/696
# gH2gr2.50year.1 <- sum(test.50year.1["H2",]>2)/696
# gH2gr4.50year.1 <- sum(test.50year.1["H2",]>4)/696
# gADgr99.50year.1 <- sum(test.50year.1["P",]>.99)/696
# gADgr95.50year.1 <- sum(test.50year.1["P",]>.95)/696

# table.1 <- signif(c(mL.50year.1,mH1.50year.1,mH2.50year.1,mAD.50year.1,
#            gH1gr2.50year.1*100,gH1gr4.50year.1*100,gH2gr2.50year.1*100,gH2gr4.50year.1*100,
#            gADgr95.50year.1*100,gADgr99.50year.1*100),3)
# names(table.1) <- c("Avg. n sites","m(H1)","m(H2)","med(p(AD))","% H1>2","% H1>4",
#                    "% H2>2","% H2>4","% p(AD)>0.95","% p(AD)>0.99")
# print(table.1)
table.1 <- c(11.30,2.94,1.67,0.99,60.90,27.00,38.20,8.62,65.80,44.70)
names(table.1) <- c("Avg. n sites","m(H1)","m(H2)","med(p(AD))","% H1>2","% H1>4",
                   "% H2>2","% H2>4","% p(AD)>0.95","% p(AD)>0.99")
print(table.1)





# Formation of groups and homogeneity.
# Classification variables:
#  2) DPLBAR, lnSAAR, lnRMED2D (weighted with regression coefficients)\n

coeff <- lm(lcv ~ DPLBAR + lnSAAR + lnRMED2D,data=potCV)$coefficients[-1]
print(coeff)

# term1 <- cd696[,"dplbar"]*coeff[1]
# term2 <- log(cd696[,"saar"])*coeff[2]
# term3 <- log(cd696[,"rmed_2d"])*coeff[3]

# roi.cd <- data.frame(cbind(term1,term2,term3))
# row.names(roi.cd) <- cd696[,"number"]

# roi02.50year <- new.env()
# for(i in 1:696) {
# print(paste(i,"/ 696"))
# assign(as.character(row.names(roi.cd)[i]), roi.st.year(roi.cd[i,],as.data.frame(roi.cd),
#         row.names(roi.cd),am696[,"am"],am696[,"number"],test="HW and AD",station.year=250,Nsim=100), env=roi02.50year)
# }
# roi02.50year <- as.list(roi02.50year)

# estrai.region <- function (x) {x$region}
# estrai.test <- function (x) {x$test}
# regioni.50year.2 <- sapply(roi02.50year, estrai.region)
# test.50year.2 <- sapply(roi02.50year, estrai.test)
# mL.50year.2 <- mean(sapply(regioni.50year.2,length))
# mH1.50year.2 <- mean(test.50year.2["H1",])
# mH2.50year.2 <- mean(test.50year.2["H2",])
# mAD.50year.2 <- median(test.50year.2["P",])
# gH1gr2.50year.2 <- sum(test.50year.2["H1",]>2)/696
# gH1gr4.50year.2 <- sum(test.50year.2["H1",]>4)/696
# gH2gr2.50year.2 <- sum(test.50year.2["H2",]>2)/696
# gH2gr4.50year.2 <- sum(test.50year.2["H2",]>4)/696
# gADgr99.50year.2 <- sum(test.50year.2["P",]>.99)/696
# gADgr95.50year.2 <- sum(test.50year.2["P",]>.95)/696

# table.2 <- signif(c(mL.50year.2,mH1.50year.2,mH2.50year.2,mAD.50year.2,
#            gH1gr2.50year.2*100,gH1gr4.50year.2*100,gH2gr2.50year.2*100,gH2gr4.50year.2*100,
#            gADgr95.50year.2*100,gADgr99.50year.2*100),3)
# names(table.2) <- c("Avg. n sites","m(H1)","m(H2)","med(p(AD))","% H1>2","% H1>4",
#                     "% H2>2","% H2>4","% p(AD)>0.95","% p(AD)>0.99")
# print(table.2)

table.2 <- c(11.20,2.95,1.64,0.99,57.20,28.60,37.90,9.34,66.40,45.10) # weigths doesnt change the result
names(table.2) <- c("Avg. n sites","m(H1)","m(H2)","med(p(AD))","% H1>2","% H1>4",
                   "% H2>2","% H2>4","% p(AD)>0.95","% p(AD)>0.99")
print(table.2)


# Formation of groups and homogeneity.
# Classification variables:
#  3) FARL, RMED2D, lnSAAR (same weight)

# term1 <- cd696[,"farl"]/sd(cd696[,"farl"])
# term2 <- cd696[,"rmed_2d"]/sd(cd696[,"rmed_2d"])
# term3 <- log(cd696[,"saar"])/sd(log(cd696[,"saar"]))

# roi.cd <- data.frame(cbind(term1,term2,term3))
# row.names(roi.cd) <- cd696[,"number"]

# roi03.50year <- new.env()
# for(i in 1:696) {
# print(paste(i,"/ 696"))
# assign(as.character(row.names(roi.cd)[i]), roi.st.year(roi.cd[i,],as.data.frame(roi.cd),
#         row.names(roi.cd),am696[,"am"],am696[,"number"],test="HW and AD",station.year=250,Nsim=100), env=roi03.50year)
# }
# roi03.50year <- as.list(roi03.50year)

# estrai.region <- function (x) {x$region}
# estrai.test <- function (x) {x$test}
# regioni.50year.3 <- sapply(roi03.50year, estrai.region)
# test.50year.3 <- sapply(roi03.50year, estrai.test)
# mL.50year.3 <- mean(sapply(regioni.50year.3,length))
# mH1.50year.3 <- mean(test.50year.3["H1",])
# mH2.50year.3 <- mean(test.50year.3["H2",])
# mAD.50year.3 <- median(test.50year.3["P",])
# gH1gr2.50year.3 <- sum(test.50year.3["H1",]>2)/696
# gH1gr4.50year.3 <- sum(test.50year.3["H1",]>4)/696
# gH2gr2.50year.3 <- sum(test.50year.3["H2",]>2)/696
# gH2gr4.50year.3 <- sum(test.50year.3["H2",]>4)/696
# gADgr99.50year.3 <- sum(test.50year.3["P",]>.99)/696
# gADgr95.50year.3 <- sum(test.50year.3["P",]>.95)/696

# table.3 <- signif(c(mL.50year.3,mH1.50year.3,mH2.50year.3,mAD.50year.3,
#            gH1gr2.50year.3*100,gH1gr4.50year.3*100,gH2gr2.50year.3*100,gH2gr4.50year.3*100,
#            gADgr95.50year.3*100,gADgr99.50year.3*100),3)
# names(table.3) <- c("Avg. n sites","m(H1)","m(H2)","med(p(AD))","% H1>2","% H1>4",
#                     "% H2>2","% H2>4","% p(AD)>0.95","% p(AD)>0.99")
# print(table.3)

table.3 <- c(11.00,3.37,1.60,0.99,67.20,35.10,34.60,8.05,67.40,46.30)
names(table.3) <- c("Avg. n sites","m(H1)","m(H2)","med(p(AD))","% H1>2","% H1>4",
                   "% H2>2","% H2>4","% p(AD)>0.95","% p(AD)>0.99")
print(table.3)




# Formation of groups and homogeneity.
# Classification variables:
#  0) lnAREA, lnSAAR, BFIHOST (same weight + lnAREA/sqrt(2))
#  1) DPLBAR, lnSAAR, lnRMED2D (with the same weight)
#  2) DPLBAR, lnSAAR, lnRMED2D (weighted with regression coefficients)
#  3) FARL, RMED2D, lnSAAR (with the same weight)\n

table.0123 <- rbind(table.0,table.1,table.2,table.3)
print(table.0123)


# Formation of groups and homogeneity (other cases).
# Classification variables:

# Classification variables:
#  0) lnAREA, lnSAAR, BFIHOST (same weight + lnAREA/sqrt(2))
#  1) DPLBAR, lnSAAR, lnRMED2D (with the same weight)
#  2) DPLBAR, lnSAAR, lnRMED2D (weighted with regression coefficients)
#  3) FARL, RMED2D, lnSAAR (with the same weight)
#  4) FARL, SMDBAR, DPLBAR, lnNGRY, lnSAAR, lnRMED2D (with the same weight)
#  5) RMED1D, ALTBAR, lnAREA, lnSAAR, lnRMED2D, lnALTBAR (with the same weight)\n

poolingroupsum <- function(roi.cd) {
 roi.50year <- new.env()
 for(i in 1:696) {
 print(paste(i,"/ 696"))
 assign(as.character(row.names(roi.cd)[i]), roi.st.year(roi.cd[i,],as.data.frame(roi.cd),
         row.names(roi.cd),am696[,"am"],am696[,"number"],test="HW and AD",station.year=250,Nsim=100), env=roi.50year)
 }
 roi.50year <- as.list(roi.50year)

 estrai.region <- function (x) {x$region}
 estrai.test <- function (x) {x$test}
 regioni.50year <- sapply(roi.50year, estrai.region)
 test.50year <- sapply(roi.50year, estrai.test)
 mL.50year <- mean(sapply(regioni.50year,length))
 mH1.50year <- mean(test.50year["H1",])
 mH2.50year <- mean(test.50year["H2",])
 mAD.50year <- median(test.50year["P",])
 gH1gr2.50year <- sum(test.50year["H1",]>2)/696
 gH1gr4.50year <- sum(test.50year["H1",]>4)/696
 gH2gr2.50year <- sum(test.50year["H2",]>2)/696
 gH2gr4.50year <- sum(test.50year["H2",]>4)/696
 gADgr99.50year <- sum(test.50year["P",]>.99)/696
 gADgr95.50year <- sum(test.50year["P",]>.95)/696
 
 table <- signif(c(mL.50year,mH1.50year,mH2.50year,mAD.50year,
            gH1gr2.50year*100,gH1gr4.50year*100,gH2gr2.50year*100,gH2gr4.50year*100,
            gADgr95.50year*100,gADgr99.50year*100),3)
 names(table) <- c("Avg. n sites","m(H1)","m(H2)","med(p(AD))","% H1>2","% H1>4",
                     "% H2>2","% H2>4","% p(AD)>0.95","% p(AD)>0.99")
 table
}

termFARL <- cd696[,"farl"]/sd(cd696[,"farl"])
termSMDBAR <- cd696[,"smdbar"]/sd(cd696[,"smdbar"])
termDPLBAR <- cd696[,"dplbar"]/sd(cd696[,"dplbar"])
termlnNGRY <- log(cd696[,"ihdtm_ngr_y"])/sd(log(cd696[,"ihdtm_ngr_y"]))
termlnSAAR <- log(cd696[,"saar"])/sd(log(cd696[,"saar"]))
termlnRMED2D <- log(cd696[,"rmed_2d"])/sd(log(cd696[,"rmed_2d"]))
termRMED1D <- cd696[,"rmed_1d"]/sd(cd696[,"rmed_1d"])
termALTBAR <- cd696[,"altbar"]/sd(cd696[,"altbar"])
termlnAREA <- log(cd696[,"dtm_area"])/sd(log(cd696[,"dtm_area"]))
termlnALTBAR <- log(cd696[,"altbar"])/sd(log(cd696[,"altbar"]))

# table.4 <- poolingroupsum(data.frame(cbind(termFARL,termSMDBAR,termDPLBAR,termlnNGRY,termlnSAAR,termlnRMED2D),
#            row.names=cd696[,"number"]))
# table.5 <- poolingroupsum(data.frame(cbind(termRMED1D,termALTBAR,termlnAREA,termlnSAAR,termlnRMED2D,termlnALTBAR),
#            row.names=cd696[,"number"]))

table.4 <- c(11.10,2.68,1.56,0.98,61.20,24.10,35.20,6.03,61.90,40.40)
names(table.4) <- c("Avg. n sites","m(H1)","m(H2)","med(p(AD))","% H1>2","% H1>4",
                   "% H2>2","% H2>4","% p(AD)>0.95","% p(AD)>0.99")
table.5 <- c(11.1,2.91,1.77,0.99,62.5,24.9,39.9,9.05,65.9,44.8)
names(table.4) <- c("Avg. n sites","m(H1)","m(H2)","med(p(AD))","% H1>2","% H1>4",
                   "% H2>2","% H2>4","% p(AD)>0.95","% p(AD)>0.99")

print(rbind(table.4,table.5))
print(rbind(table.0123,table.4,table.5))

# ------------------------------------------------------------------------------------- #
# ------------------------------------------------------------------------------------- #
#                                       THE END                                         #
# ------------------------------------------------------------------------------------- #
# ------------------------------------------------------------------------------------- #

