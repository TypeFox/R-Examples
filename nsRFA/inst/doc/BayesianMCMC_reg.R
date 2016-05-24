### R code from vignette source 'BayesianMCMC_reg.Rnw'

###################################################
### code chunk number 1: BayesianMCMC_reg.Rnw:51-54
###################################################
library(nsRFA)
data(Ardechedata)
ls()


###################################################
### code chunk number 2: BayesianMCMC_reg.Rnw:59-60
###################################################
str(SaintMartin_cont)


###################################################
### code chunk number 3: BayesianMCMC_reg.Rnw:64-65
###################################################
str(SaintMartin_hist)


###################################################
### code chunk number 4: plot_SaintMartin_cont_hist (eval = FALSE)
###################################################
## plot(c(1640,2010), c(0, 8000), type="n", xlab="", ylab="discharge (m3/s)")
##  points(SaintMartin_cont[,"year"], SaintMartin_cont[,"peak"], type="o", pch=21, bg="white")
##  points(SaintMartin_hist$peaks[,"year"], SaintMartin_hist$peaks[,"peak"], pch=16)
##  points(SaintMartin_hist$peaks[,"year"], SaintMartin_hist$peaks[,"peak"], type="h")
##  segments(x0=SaintMartin_hist$thresholds[,"from.yr"], 
##           x1=SaintMartin_hist$thresholds[,"to.yr"],
##           y0=SaintMartin_hist$thresholds[,"threshold"], lty=2)


###################################################
### code chunk number 5: BayesianMCMC_reg.Rnw:79-81
###################################################
par(mar=c(3,3,2,1)+0.03, mgp=c(1.5,0.3,0), tcl=.2, xaxs="i", yaxs="i", las=0)
plot(c(1640,2010), c(0, 8000), type="n", xlab="", ylab="discharge (m3/s)")
 points(SaintMartin_cont[,"year"], SaintMartin_cont[,"peak"], type="o", pch=21, bg="white")
 points(SaintMartin_hist$peaks[,"year"], SaintMartin_hist$peaks[,"peak"], pch=16)
 points(SaintMartin_hist$peaks[,"year"], SaintMartin_hist$peaks[,"peak"], type="h")
 segments(x0=SaintMartin_hist$thresholds[,"from.yr"], 
          x1=SaintMartin_hist$thresholds[,"to.yr"],
          y0=SaintMartin_hist$thresholds[,"threshold"], lty=2)


###################################################
### code chunk number 6: BayesianMCMC_reg.Rnw:92-96
###################################################
xhistSM <- c(-1,SaintMartin_hist$peaks[,"peak"])
seuilSM <- c(7250, 6000, rep(5050, 10), rep(2400, 21))
nbansSM <- c(127, 55, 65, rep(0,9), 71, rep(0, 20))
data.frame(xhistSM, seuilSM, nbansSM)


###################################################
### code chunk number 7: BayesianMCMC_reg.Rnw:105-106
###################################################
set.seed(198)


###################################################
### code chunk number 8: BayesianMCMC_reg.Rnw:108-114
###################################################
fit <- BayesianMCMC(xcont=SaintMartin_cont[,"peak"], 
                    xhist=xhistSM, 
                    seuil=seuilSM, 
                    nbans=nbansSM, 
                    nbpas=1000, nbchaines=3, confint=c(0.05, 0.95), dist="GEV",
                    varparameters0=c(NA, NA, 0.5))


###################################################
### code chunk number 9: BayesianMCMC_reg.Rnw:120-121
###################################################
fit <- BayesianMCMCcont(fit)


###################################################
### code chunk number 10: BayesianMCMC_reg.Rnw:125-126 (eval = FALSE)
###################################################
## plot(fit, which=1:4, ask=TRUE)


###################################################
### code chunk number 11: plot_fit_1
###################################################
par(mar=c(3,3,2,1)+0.03, mgp=c(1.5,0.3,0), tcl=.2, xaxs="r", yaxs="r", las=0)
plot(fit, 1)


###################################################
### code chunk number 12: BayesianMCMC_reg.Rnw:145-147
###################################################
par(mar=c(3,3,2,1)+0.03, mgp=c(1.5,0.3,0), tcl=.2, xaxs="r", yaxs="r", las=0)
plot(fit, 2)


###################################################
### code chunk number 13: plot_fit_3
###################################################
par(mar=c(3,3,2,1)+0.03, mgp=c(1.5,0.3,0), tcl=.2, xaxs="r", yaxs="r", las=0)
plot(fit, 3)


###################################################
### code chunk number 14: plot_fit_4
###################################################
par(mar=c(3,3,2,1)+0.03, mgp=c(1.5,0.3,0), tcl=.2, xaxs="r", yaxs="r", las=0)
plot(fit, 4)


###################################################
### code chunk number 15: BayesianMCMC_reg.Rnw:184-189 (eval = FALSE)
###################################################
## SaintMartin_cont
## Vogue_cont
## SaintLaurent_cont
## Beauvene_cont
## Chambonas_cont


###################################################
### code chunk number 16: plot_Ardeche_cont (eval = FALSE)
###################################################
## plot(SaintMartin_cont[,"year"], SaintMartin_cont[,"peak"], type="o", pch=21, bg="white", 
##      xlab="", ylab="discharge (m3/s)", main=" Saint Martin on the Ardeche River (2240 km2)")
## plot(Vogue_cont[,"year"], Vogue_cont[,"peak"], type="o", pch=21, bg="white", 
##      xlab="", ylab="discharge (m3/s)", main="Vogue on the Ardeche River (636 km2)")
## plot(SaintLaurent_cont[,"year"], SaintLaurent_cont[,"peak"], type="o", pch=21, bg="white", 
##      xlab="", ylab="discharge (m3/s)", main="Saint Laurent on the Borne River (63 km2)")
## plot(Beauvene_cont[,"year"], Beauvene_cont[,"peak"], type="o", pch=21, bg="white", 
##      xlab="", ylab="discharge (m3/s)", main="Beauvene on the Eyrieux River (392 km2)")
## plot(Chambonas_cont[,"year"], Chambonas_cont[,"peak"], type="o", pch=21, bg="white", 
##      xlab="", ylab="discharge (m3/s)", main="Chambonas on the Chassezac River (507 km2)")


###################################################
### code chunk number 17: BayesianMCMC_reg.Rnw:204-207
###################################################
layout(matrix(1:6, ncol=2))
par(mar=c(3,3,2,1)+0.03, mgp=c(1.5,0.3,0), tcl=.2, xaxs="r", yaxs="r", las=0)
plot(SaintMartin_cont[,"year"], SaintMartin_cont[,"peak"], type="o", pch=21, bg="white", 
     xlab="", ylab="discharge (m3/s)", main=" Saint Martin on the Ardeche River (2240 km2)")
plot(Vogue_cont[,"year"], Vogue_cont[,"peak"], type="o", pch=21, bg="white", 
     xlab="", ylab="discharge (m3/s)", main="Vogue on the Ardeche River (636 km2)")
plot(SaintLaurent_cont[,"year"], SaintLaurent_cont[,"peak"], type="o", pch=21, bg="white", 
     xlab="", ylab="discharge (m3/s)", main="Saint Laurent on the Borne River (63 km2)")
plot(Beauvene_cont[,"year"], Beauvene_cont[,"peak"], type="o", pch=21, bg="white", 
     xlab="", ylab="discharge (m3/s)", main="Beauvene on the Eyrieux River (392 km2)")
plot(Chambonas_cont[,"year"], Chambonas_cont[,"peak"], type="o", pch=21, bg="white", 
     xlab="", ylab="discharge (m3/s)", main="Chambonas on the Chassezac River (507 km2)")


###################################################
### code chunk number 18: BayesianMCMC_reg.Rnw:211-213
###################################################
Ardeche_areas  # km2
Ardeche_ungauged_extremes  # m3/s


###################################################
### code chunk number 19: BayesianMCMC_reg.Rnw:216-217
###################################################
Ardeche_ungauged_extremes_new <- merge(Ardeche_ungauged_extremes, Ardeche_areas)


###################################################
### code chunk number 20: BayesianMCMC_reg.Rnw:232-246
###################################################
xcont <- c(SaintMartin_cont[,2],
           Vogue_cont[,2], 
           SaintLaurent_cont[,2], 
           Beauvene_cont[,2], Chambonas_cont[,2])
scont <- c(rep(2240, length(SaintMartin_cont[,2])), 
           rep(636, length(Vogue_cont[,2])), 
           rep(63, length(SaintLaurent_cont[,2])), 
           rep(392, length(Beauvene_cont[,2])), 
           rep(507, length(Chambonas_cont[,2])))                           
xhist <- Ardeche_ungauged_extremes_new[,"peak"]
xhist[Ardeche_ungauged_extremes_new[,"station"]=="Chambonas"] <- -1
shist <- Ardeche_ungauged_extremes_new[,"area"]
nbans <- Ardeche_ungauged_extremes_new[,"associated.period"]
seuil <- Ardeche_ungauged_extremes_new[,"peak"]


###################################################
### code chunk number 21: BayesianMCMC_reg.Rnw:250-251
###################################################
set.seed(198)


###################################################
### code chunk number 22: BayesianMCMC_reg.Rnw:253-263
###################################################
fitreg <- BayesianMCMCreg(xcont=xcont, 
                          scont=scont,
                          xhist=xhist,
                          shist=shist,
                          nbans=nbans,
                          seuil=seuil,
                          nbpas=1000, nbchaines=3, confint=c(0.05,0.95), dist="GEV",
                          varparameters0=c(NA, NA, NA, 0.5))
fitreg <- BayesianMCMCregcont(fitreg)
fitreg <- BayesianMCMCregcont(fitreg)


###################################################
### code chunk number 23: BayesianMCMC_reg.Rnw:281-282 (eval = FALSE)
###################################################
## plot(fitreg, which=1:4, ask=TRUE)


###################################################
### code chunk number 24: plot_fitreg_1
###################################################
par(mar=c(3,3,2,1)+0.03, mgp=c(1.5,0.3,0), tcl=.2, xaxs="r", yaxs="r", las=0)
plot(fitreg, 1)


###################################################
### code chunk number 25: BayesianMCMC_reg.Rnw:295-297
###################################################
par(mar=c(3,3,2,1)+0.03, mgp=c(1.5,0.3,0), tcl=.2, xaxs="r", yaxs="r", las=0)
plot(fitreg, 2)


###################################################
### code chunk number 26: plot_fitreg_3
###################################################
par(mar=c(3,3,2,1)+0.03, mgp=c(1.5,0.3,0), tcl=.2, xaxs="r", yaxs="r", las=0)
plot(fitreg, 3)


###################################################
### code chunk number 27: plot_fitreg_4
###################################################
par(mar=c(3,3,2,1)+0.03, mgp=c(1.5,0.3,0), tcl=.2, xaxs="r", yaxs="r", las=0)
plot(fitreg, 4)


###################################################
### code chunk number 28: plotBayesianMCMCreg_surfArdeche (eval = FALSE)
###################################################
## plotBayesianMCMCreg_surf(fitreg, surf=unique(Ardeche_areas$area))


###################################################
### code chunk number 29: plotBayesianMCMCreg_surfArdechegrosso
###################################################
layout(matrix(1:20, ncol=4, byrow=TRUE))
par(mar=c(3,3,2,1)+0.03, mgp=c(1.5,0.3,0), tcl=.2, xaxs="r", yaxs="r", las=0)
plotBayesianMCMCreg_surf(fitreg, surf=unique(Ardeche_areas$area))


