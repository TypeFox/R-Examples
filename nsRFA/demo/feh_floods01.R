# ------------------------------------------------------------------------------------- #
# ------------------------------------------------------------------------------------- #
#                        ANALYSIS OF FEH DATA: INDEX-VALUE                              #
# ------------------------------------------------------------------------------------- #
# ------------------------------------------------------------------------------------- #


# Data loading

data(FEH1000)


# Sites used in model development (pag.102 FEH Vol.3):
#  area>0.5 km2
#  digital catchment data available
#  urbext<0.025\n

urbext <- cd[,"urbext1990"]
area <- cd[,"dtm_area"]
cd732 <- cd[(!is.nan(cd[,"dtm_area"]))&(urbext<0.025)&(area>0.5),] # vs 687 - 728 of FEH

fac <- factor(am[,"number"],levels=cd732[,"number"])
am732 <- am[!is.na(fac),]
nlevels(as.factor(am732[,"number"]))



# Index-flood = median

QMED <- tapply(am732[,4],am732[,1],median)
lnQMED <- log(QMED)

# Catchment descriptors (fig. 13.1 pag. 104)

lnAREA <- log(cd732[,"dtm_area"])
lnDPLBAR <- log(cd732[,"dplbar"])
lnSPRHOST <- log(cd732[,"sprhost"])
lnBFIHOST <- log(cd732[,"bfihost"])
lnSAAR <- log(cd732[,"saar"])
lnRMED1 <- log(cd732[,"rmed_1d"])
#lnNWET <- log(cd732[,""])
lnALTBAR <- log(cd732[,"altbar"])
lnFARL <- log(cd732[,"farl"])

M <- data.frame(cbind(lnQMED,lnAREA,lnDPLBAR,lnSPRHOST,lnBFIHOST,lnSAAR,lnRMED1,lnALTBAR,lnFARL))
print(cor(M))
plot(M,pch=".",cex=2)
par(ask = interactive())

# Additional variables (pag. 105):

RESHOST <- cd732[,"bfihost"] + 1.30*(cd732[,"sprhost"]/100)-0.987
lnAREAsq <- lnAREA^2
lnSAARsq <- lnSAAR^2

M <- data.frame(cbind(M,RESHOST,lnAREAsq,lnSAARsq))



# Ordinary Least Square models:
# create a function using the 'leaps' function of package 'subselect'
# to perform all-possible-regressions: 
#   bestregressions <- function(dip,ind) {
#    Y <- as.numeric(dip)
#    X <- ind
#    Sy <- var(Y)
#    Sx <- var(X)
#    Sxy <- var(X,Y)
#    Dm.mat <- Sx
#    Dm.H <- Sxy %*% t(Sxy)/Sy
#    require(subselect)
#    Dm.leaps <- leaps(Dm.mat, kmin=1, kmax=6, H=Dm.H, r=1, nsol=3)
#    Dm.leaps
#    for(i in 6:1) {for(j in 1:3) {print(colnames(X)[Dm.leaps$subsets[j,c(1:6),i]])}}
#   }

#   bestregressions(M[,1],M[,-1])
# [1] "lnAREA"    "lnSPRHOST" "lnSAAR"    "lnFARL"    "RESHOST"   "lnSAARsq"
# [1] "lnAREA"    "lnSPRHOST" "lnBFIHOST" "lnSAAR"    "lnFARL"    "lnSAARsq"
# [1] "lnAREA"    "lnSPRHOST" "lnSAAR"    "lnRMED1"   "lnFARL"    "lnSAARsq"
# [1] "lnAREA"    "lnSPRHOST" "lnSAAR"    "lnFARL"    "lnSAARsq"
# [1] "lnAREA"    "lnSPRHOST" "lnSAAR"    "lnFARL"    "RESHOST"
# [1] "lnAREA"    "lnSPRHOST" "lnSAAR"    "RESHOST"   "lnSAARsq"
# [1] "lnAREA"    "lnSPRHOST" "lnSAAR"    "lnSAARsq"
# [1] "lnAREA"    "lnSPRHOST" "lnSAAR"    "lnFARL"
# [1] "lnAREA"    "lnSPRHOST" "lnSAAR"    "RESHOST"
# [1] "lnAREA"    "lnSPRHOST" "lnSAAR"
# [1] "lnAREA"    "lnSPRHOST" "lnSAARsq"
# [1] "lnAREA"    "lnBFIHOST" "lnSAAR"
# [1] "lnAREA" "lnSAAR"
# [1] "lnAREA"   "lnSAARsq"
# [1] "lnAREA"    "lnSPRHOST"
# [1] "lnAREAsq"
# [1] "lnAREA"
# [1] "lnDPLBAR"


# Ordinary Least Square models (graphics and statistics):

graphics.lm <- function (regr) {
 par(mfrow=c(2,2), cex=.7)
 plot(regr$fitted.values,regr$residuals,pch=".",cex=3,xlab="lnQMED Fitted",ylab="lnQMED Residuals")
 abline(0,0,lty=3)
 normplot(regr$residuals,pch=".",cex=3,xlab="lnQMED Residuals")
 
 plot(regr$fitted.values,lnQMED,pch=".",cex=3,xlab="lnQMED Originals",ylab="lnQMED Fitted")
 abline(0,1,lty=3)
 intervals <- predinterval.lm(regr)
 intervals <- intervals[order(intervals[,1]),]
 lines(intervals[,c(1,2)],lty=2)
 lines(intervals[,c(1,3)],lty=2)
 Rsq <- signif(R2.lm(regr),3)
 rmse <- signif(RMSE.lm(regr),3)
 rmsep <- signif(RMSEP(lnQMED,regr$fitted.values)*100,3)
 mtext(paste("R2 = ",Rsq),3,-1.5,adj=0.05,cex=.7)
 mtext(paste("RMSE = ",rmse),1,-2.5,adj=0.95,cex=.7)
 mtext(paste("RMSEP = ",rmsep,"%"),1,-1.5,adj=0.95,cex=.7)

 plot(exp(regr$fitted.values),exp(lnQMED),pch=".",cex=3,xlab="QMED Originals",ylab="QMED Fitted")
 abline(0,1,lty=3)
 lines(exp(intervals[,c(1,2)]),lty=2)
 lines(exp(intervals[,c(1,3)]),lty=2)
 Rsq <- signif(R2(exp(lnQMED),exp(regr$fitted.values)),3)
 rmse <- signif(RMSE(exp(lnQMED),exp(regr$fitted.values)),3)
 rmsep <- signif(RMSEP(exp(lnQMED),exp(regr$fitted.values))*100,3)
 mtext(paste("R2 = ",Rsq),3,-1.5,adj=0.05,cex=.7)
 mtext(paste("RMSE = ",rmse),1,-2.5,adj=0.95,cex=.7)
 mtext(paste("RMSEP = ",rmsep,"%"),1,-1.5,adj=0.95,cex=.7)

 par(mfrow=c(1,1),cex=1)
 title(main=paste(names(regr$coefficients)[-1], collapse=", "),cex.main=.7,font.main=1)
}

graphics.lm(lm(lnQMED ~ lnDPLBAR))
graphics.lm(lm(lnQMED ~ lnAREA))
graphics.lm(lm(lnQMED ~ lnAREAsq))
graphics.lm(lm(lnQMED ~ lnAREA + lnSPRHOST))
graphics.lm(lm(lnQMED ~ lnAREA + lnSAARsq))
graphics.lm(lm(lnQMED ~ lnAREA + lnSAAR))
graphics.lm(lm(lnQMED ~ lnAREA + lnBFIHOST + lnSAAR))
graphics.lm(lm(lnQMED ~ lnAREA + lnSPRHOST + lnSAARsq))
graphics.lm(lm(lnQMED ~ lnAREA + lnSPRHOST + lnSAAR))
graphics.lm(lm(lnQMED ~ lnAREA + lnSPRHOST + lnSAAR + RESHOST))
graphics.lm(lm(lnQMED ~ lnAREA + lnSPRHOST + lnSAAR + lnFARL))
graphics.lm(lm(lnQMED ~ lnAREA + lnSPRHOST + lnSAAR + lnSAARsq))
graphics.lm(lm(lnQMED ~ lnAREA + lnSPRHOST + lnSAAR + RESHOST + lnSAARsq))
graphics.lm(lm(lnQMED ~ lnAREA + lnSPRHOST + lnSAAR + lnFARL + RESHOST))
graphics.lm(lm(lnQMED ~ lnAREA + lnSPRHOST + lnSAAR + lnFARL + lnSAARsq))
graphics.lm(lm(lnQMED ~ lnAREA + lnSPRHOST + lnSAAR + lnRMED1 + lnFARL + lnSAARsq))
graphics.lm(lm(lnQMED ~ lnAREA + lnSPRHOST + lnBFIHOST + lnSAAR + lnFARL + lnSAARsq))
graphics.lm(lm(lnQMED ~ lnAREA + lnSPRHOST + lnSAAR + lnFARL + RESHOST + lnSAARsq))



# ------------------------------------------------------------------------------------- #
# ------------------------------------------------------------------------------------- #
#                                       THE END                                         # 
# ------------------------------------------------------------------------------------- #
# ------------------------------------------------------------------------------------- #

