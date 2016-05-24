### R code from vignette source 'ex-npbr.rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: ex-npbr.rnw:60-65
###################################################
owidth <- getOption("width")
options("width"=70)
ow <- getOption("warn")
options("warn"=-1)
.PngNo <- 0


###################################################
### code chunk number 2: bfig (eval = FALSE)
###################################################
## .PngNo <- .PngNo + 1; name.file <- paste("Fig-bitmap-", .PngNo, ".pdf", sep="")
## pdf(file=name.file, width = 9, height = 7, pointsize = 14, bg = "white")


###################################################
### code chunk number 3: bfig2 (eval = FALSE)
###################################################
## .PngNo <- .PngNo + 1; name.file <- paste("Fig-bitmap-", .PngNo, ".pdf", sep="")
## pdf(file=name.file, width = 14, height = 7, pointsize = 14, bg = "white")


###################################################
### code chunk number 4: bfig3 (eval = FALSE)
###################################################
## .PngNo <- .PngNo + 1; name.file <- paste("Fig-bitmap-", .PngNo, ".pdf", sep="")
## pdf(file=name.file, width = 7, height = 7, pointsize = 14, bg = "white")


###################################################
### code chunk number 5: bfig4 (eval = FALSE)
###################################################
## .PngNo <- .PngNo + 1; name.file <- paste("Fig-bitmap-", .PngNo, ".pdf", sep="")
## pdf(file=name.file, width = 18, height = 7, pointsize = 14, bg = "white")


###################################################
### code chunk number 6: bfig5 (eval = FALSE)
###################################################
## .PngNo <- .PngNo + 1; name.file <- paste("Fig-bitmap-", .PngNo, ".pdf", sep="")
## pdf(file=name.file, width = 18, height = 14, pointsize = 14, bg = "white")


###################################################
### code chunk number 7: zfig (eval = FALSE)
###################################################
## dev.null <- dev.off()
## cat("\\includegraphics[width=0.75\\textwidth]{", name.file, "}\n\n", sep="")


###################################################
### code chunk number 8: zfig2 (eval = FALSE)
###################################################
## dev.null <- dev.off()
## cat("\\includegraphics[width=0.9\\textwidth]{", name.file, "}\n\n", sep="")


###################################################
### code chunk number 9: zfig3 (eval = FALSE)
###################################################
## dev.null <- dev.off()
## cat("\\includegraphics[width=0.8\\textwidth]{", name.file, "}\n\n", sep="")


###################################################
### code chunk number 10: ex-npbr.rnw:357-363 (eval = FALSE)
###################################################
## require("npbr")
## data("records")
## data("nuclear")
## data("air")
## data("post")
## data("green")


###################################################
### code chunk number 11: ex-npbr.rnw:366-375 (eval = FALSE)
###################################################
## plot(result~year, data=records, col='blue2', 
##  xlab="year", ylab="1500m record")
## plot(ytab~xtab, data=nuclear, col='blue2', 
##  xlab="temp. of the reactor vessel", ylab="fracture toughness")
## plot(ytab~xtab, data=air, col='blue2', 
##  xlab="input", ylab="output") 
## plot(yprod~xinput, data=post, col='blue2', 
## xlab="quantity of labor", ylab="volume of delivered mail")
## plot(log(OUTPUT)~log(COST), data=green, col='blue2')


###################################################
### code chunk number 12: ex-npbr.rnw:380-394 (eval = FALSE)
###################################################
## .PngNo <- .PngNo + 1; name.file <- paste("Fig-bitmap-", .PngNo, ".pdf", sep="")
## pdf(file=name.file, width = 18, height = 14, pointsize = 14, bg = "white")
## op<-par(mfrow=c(2,3),mar=c(3,3.1,2.1,2.1),mgp=c(2,.4,0),oma=c(0,0,0,0),
##  cex.lab=1.2, cex.main=1, col.main="blue")
## plot(result~year, data=records, col='blue2', pch=1,
##  xlab="year", ylab="1500m record", main="(a)") 
## plot(ytab~xtab, data=nuclear, pch=1,col='blue2', 
##  xlab="temp. of the reactor vessel", ylab="fracture toughness", main="(b)") 
## plot(ytab~xtab, data=air, pch=1,col='blue2', 
##  xlab="input", ylab="output", main="(c)")  
## plot(yprod~xinput, data=post, pch=1, col='blue2', 
## xlab="quantity of labor", ylab="volume of delivered mail", main="(d)") 
## plot(log(OUTPUT)~log(COST), data=green, pch=1,col='blue2', main="(e)") 
## par(op)
## dev.null <- dev.off()
## cat("\\includegraphics[width=0.9\\textwidth]{", name.file, "}\n\n", sep="")


###################################################
### code chunk number 13: ex-npbr.rnw:465-468 (eval = FALSE)
###################################################
## x.air <- seq(min(air$xtab), max(air$xtab), length.out=101)
## x.green <- seq(min(log(green$COST)), max(log(green$COST)), 
##  length.out=101)


###################################################
### code chunk number 14: ex-npbr.rnw:471-481 (eval = FALSE)
###################################################
## y.dea.green<-dea_est(log(green$COST), log(green$OUTPUT), 
##  x.green, type="dea")
## y.fdh.green<-dea_est(log(green$COST), log(green$OUTPUT), 
##  x.green, type="fdh")
## y.lfdh.green=dea_est(log(green$COST), log(green$OUTPUT), 
##  x.green, type="lfdh")
##  
## y.dea.air<-dea_est(air$xtab, air$ytab, x.air, type="dea")
## y.fdh.air<-dea_est(air$xtab, air$ytab, x.air, type="fdh")
## y.lfdh.air=dea_est(air$xtab, air$ytab, x.air, type="lfdh") 


###################################################
### code chunk number 15: ex-npbr.rnw:485-500 (eval = FALSE)
###################################################
## plot(x.green, y.dea.green,  lty=4, lwd=4, col="cyan",
##  type="l", xlab="log(cost)",ylab="log(output)")    
## lines(x.green, y.fdh.green, lty=1, lwd=4, col="green")
## lines(x.green, y.lfdh.green, lty=2, lwd=4, col="magenta")   
## legend("topleft", legend=c("DEA","FDH","LFDH"), 
##  col=c("cyan","green","magenta"), lty=c(4,1,2), lwd=4)
## points(log(OUTPUT)~log(COST), data=green, cex=1) 
##  
## plot(x.air, y.dea.air,  lty=4, lwd=4, col="cyan",
##  type="l", xlab="input",ylab="output")     
## lines(x.air, y.fdh.air, lty=1, lwd=4, col="green")
## lines(x.air, y.lfdh.air, lty=2, lwd=4, col="magenta")   
## legend("topleft", legend=c("DEA","FDH","LFDH"), 
##  col=c("cyan","green","magenta"), lty=c(4,1,2), lwd=4)
## points(ytab~xtab, data=air)  


###################################################
### code chunk number 16: ex-npbr.rnw:504-523 (eval = FALSE)
###################################################
## .PngNo <- .PngNo + 1; name.file <- paste("Fig-bitmap-", .PngNo, ".pdf", sep="")
## pdf(file=name.file, width = 14, height = 7, pointsize = 14, bg = "white")
## op=par(mfrow=c(1,2),mar=c(3,3.1,2.1,2.1),mgp=c(2,.4,0),oma=c(0,0,0,0),cex.lab=1.2)
## plot(x.green, y.dea.green,  lty=4, lwd=4, col="cyan",
##  type="l", xlab="log(cost)",ylab="log(output)")    
## lines(x.green, y.fdh.green, lty=1, lwd=4, col="green")
## lines(x.green, y.lfdh.green, lty=2, lwd=4, col="magenta")   
## legend("topleft", legend=c("DEA","FDH","LFDH"), 
##  col=c("cyan","green","magenta"), lty=c(4,1,2), lwd=4)
## points(log(OUTPUT)~log(COST), data=green, cex=1) 
##  
## plot(x.air, y.dea.air,  lty=4, lwd=4, col="cyan",
##  type="l", xlab="input",ylab="output")     
## lines(x.air, y.fdh.air, lty=1, lwd=4, col="green")
## lines(x.air, y.lfdh.air, lty=2, lwd=4, col="magenta")   
## legend("topleft", legend=c("DEA","FDH","LFDH"), 
##  col=c("cyan","green","magenta"), lty=c(4,1,2), lwd=4)
## points(ytab~xtab, data=air)  
## par(op) 
## dev.null <- dev.off()
## cat("\\includegraphics[width=0.9\\textwidth]{", name.file, "}\n\n", sep="")


###################################################
### code chunk number 17: ex-npbr.rnw:562-567 (eval = FALSE)
###################################################
## (p.aic.records<-poly_degree(records$year, 1/records$result, prange=0:12, 
##  type = "AIC"))
## (p.aic.air<-poly_degree(air$xtab, air$ytab, 
##  type = "AIC"))
## (p.aic.nuc<-poly_degree(nuclear$xtab, nuclear$ytab, type = "AIC"))


###################################################
### code chunk number 18: ex-npbr.rnw:571-580 (eval = FALSE)
###################################################
## x.records<-seq(min(records$year), max(records$year), length.out=101)
## y.poly.records<-poly_est(records$year, 1/records$result, x.records, 
##  deg=p.aic.records)
## y.poly.air<-poly_est(air$xtab, air$ytab, x.air, 
##  deg=p.aic.air)
## x.nucl <- seq(min(nuclear$xtab), max(nuclear$xtab), 
##  length.out=101) 
## y.poly.nuc<-poly_est(nuclear$xtab, nuclear$ytab, x.nucl, 
##  deg=p.aic.nuc) 


###################################################
### code chunk number 19: ex-npbr.rnw:586-599 (eval = FALSE)
###################################################
## plot(x.records, 1/y.poly.records, lty=1, lwd=4, 
##  col="magenta", type="l")
## points(result~year, data=records) 
## plot(x.air, y.poly.air, lty=1, lwd=4, 
##  col="magenta", type="l")
## points(ytab~xtab, data=air)  
## legend("topleft",legend=paste("degree =",p.aic.air), 
##  col="magenta", lwd=4, lty=1) 
## plot(x.nucl, y.poly.nuc, lty=1, lwd=4, 
##  col="cyan", type="l", ylim=range(nuclear$ytab))
## points(ytab~xtab, data=nuclear)  
## legend("topleft",legend=paste("degree =",p.aic.nuc), 
##  col="cyan", lwd=4, lty=1)


###################################################
### code chunk number 20: ex-npbr.rnw:603-622 (eval = FALSE)
###################################################
## .PngNo <- .PngNo + 1; name.file <- paste("Fig-bitmap-", .PngNo, ".pdf", sep="")
## pdf(file=name.file, width = 18, height = 7, pointsize = 14, bg = "white")
## op=par(mfrow=c(1,3),mar=c(3,3.1,2.1,2.1),mgp=c(2,.4,0),oma=c(0,0,0,0),cex.lab=1.2)
## plot(x.records, 1/y.poly.records, lty=1, lwd=4, 
##  col="green", type="l",  xlab="year", ylab="1500m record")
## points(result~year, data=records) 
## legend("topleft",legend=paste("degree =",p.aic.records), 
##  col="green", lwd=4, lty=1)
## plot(x.air, y.poly.air, lty=1, lwd=4, 
##  col="magenta", type="l", xlab="input", ylab="output")
## points(ytab~xtab, data=air)  
## legend("topleft",legend=paste("degree =",p.aic.air), 
##  col="magenta", lwd=4, lty=1) 
## plot(x.nucl, y.poly.nuc, lty=1, lwd=4, 
##  col="cyan", type="l", ylim=range(nuclear$ytab), xlab="temperature", ylab="toughness")
## points(ytab~xtab, data=nuclear)  
## legend("topleft",legend=paste("degree =",p.aic.nuc), 
##  col="cyan", lwd=4, lty=1)
## par(op)
## dev.null <- dev.off()
## cat("\\includegraphics[width=0.9\\textwidth]{", name.file, "}\n\n", sep="")


###################################################
### code chunk number 21: ex-npbr.rnw:726-730 (eval = FALSE)
###################################################
## (kn.bic.air.u<-quad_spline_kn(air$xtab, 
##  air$ytab, method="u", type="BIC"))
## (kn.bic.green.u<-quad_spline_kn(log(green$COST), 
##  log(green$OUTPUT), method="u", type="BIC"))


###################################################
### code chunk number 22: ex-npbr.rnw:734-738 (eval = FALSE)
###################################################
## y.quad.air.u<-quad_spline_est(air$xtab, 
##  air$ytab, x.air, kn=kn.bic.air.u, method="u")
## y.quad.green.u<-quad_spline_est(log(green$COST), 
##  log(green$OUTPUT), x.green, kn=kn.bic.green.u, method="u")


###################################################
### code chunk number 23: ex-npbr.rnw:741-745 (eval = FALSE)
###################################################
## (kn.bic.air.m<-quad_spline_kn(air$xtab, 
##  air$ytab, method="m", type="BIC"))
## (kn.bic.green.m<-quad_spline_kn(log(green$COST), 
##  log(green$OUTPUT), method="m", type="BIC"))


###################################################
### code chunk number 24: ex-npbr.rnw:749-753 (eval = FALSE)
###################################################
## y.quad.air.m<-quad_spline_est(air$xtab, 
##  air$ytab, x.air, kn=kn.bic.air.m, method="m")
## y.quad.green.m<-quad_spline_est(log(green$COST), 
##  log(green$OUTPUT), x.green, kn=kn.bic.green.m, method="m")


###################################################
### code chunk number 25: ex-npbr.rnw:758-762 (eval = FALSE)
###################################################
## (kn.bic.air.mc<-quad_spline_kn(air$xtab, 
##  air$ytab, method="mc", type="BIC"))
## (kn.bic.green.mc<-quad_spline_kn(log(green$COST), 
##  log(green$OUTPUT), method="mc", type="BIC"))


###################################################
### code chunk number 26: ex-npbr.rnw:766-771 (eval = FALSE)
###################################################
## y.quad.air.mc<-quad_spline_est(air$xtab, air$ytab, x.air, 
##  kn=kn.bic.air.mc, method="mc", all.dea=TRUE)
## y.quad.green.mc<-quad_spline_est(log(green$COST), 
##  log(green$OUTPUT), x.green, kn=kn.bic.green.mc, 
##  method="mc", all.dea=TRUE)


###################################################
### code chunk number 27: ex-npbr.rnw:776-792 (eval = FALSE)
###################################################
## plot(x.air, y.quad.air.u, lty=1, lwd=4, col="green", 
##  type="l", xlab="input", ylab="output")   
## lines(x.air, y.quad.air.m, lty=2, lwd=4, col="cyan")
## lines(x.air, y.quad.air.mc, lty=3, lwd=4, col="magenta")   
## points(ytab~xtab, data=air)
## legend("topleft", col=c("green","cyan","magenta"), 
##  lty=c(1,2,3), legend=c("unconstrained", "monotone", 
##  "monotone + concave"), lwd=4, cex=0.8)
## plot(x.green, y.quad.green.u, lty=1, lwd=4, col="green", 
##  type="l", xlab="log(COST)", ylab="log(OUTPUT)")   
## lines(x.green, y.quad.green.m, lty=2, lwd=4, col="cyan")
## lines(x.green, y.quad.green.mc, lwd=4, lty=3, col="magenta")   
## points(log(OUTPUT)~log(COST), data=green)
## legend("topleft", col=c("green","cyan","magenta"), 
## lty=c(1,2,3), legend=c("unconstrained", "monotone", 
##  "monotone + concave"), lwd=4, cex=0.8)


###################################################
### code chunk number 28: ex-npbr.rnw:796-816 (eval = FALSE)
###################################################
## .PngNo <- .PngNo + 1; name.file <- paste("Fig-bitmap-", .PngNo, ".pdf", sep="")
## pdf(file=name.file, width = 14, height = 7, pointsize = 14, bg = "white")
## op=par(mfrow=c(1,2),mar=c(3,3.1,2.1,2.1),mgp=c(2,.4,0),oma=c(0,0,0,0),cex.lab=1.2)
## plot(x.air, y.quad.air.u, lty=1, lwd=4, col="green", 
##  type="l", xlab="input", ylab="output")   
## lines(x.air, y.quad.air.m, lty=2, lwd=4, col="cyan")
## lines(x.air, y.quad.air.mc, lty=3, lwd=4, col="magenta")   
## points(ytab~xtab, data=air)
## legend("topleft", col=c("green","cyan","magenta"), 
##  lty=c(1,2,3), legend=c("unconstrained", "monotone", 
##  "monotone + concave"), lwd=4, cex=0.8)
## plot(x.green, y.quad.green.u, lty=1, lwd=4, col="green", 
##  type="l", xlab="log(COST)", ylab="log(OUTPUT)")   
## lines(x.green, y.quad.green.m, lty=2, lwd=4, col="cyan")
## lines(x.green, y.quad.green.mc, lwd=4, lty=3, col="magenta")   
## points(log(OUTPUT)~log(COST), data=green)
## legend("topleft", col=c("green","cyan","magenta"), 
## lty=c(1,2,3), legend=c("unconstrained", "monotone", 
##  "monotone + concave"), lwd=4, cex=0.8)
## par(op)
## dev.null <- dev.off()
## cat("\\includegraphics[width=0.9\\textwidth]{", name.file, "}\n\n", sep="")


###################################################
### code chunk number 29: ex-npbr.rnw:851-863 (eval = FALSE)
###################################################
## (kn.bic.air.u<-cub_spline_kn(air$xtab, air$ytab, 
##   method="u", type="BIC"))
## (kn.bic.green.u<-cub_spline_kn(log(green$COST), 
##   log(green$OUTPUT), method="u", type="BIC"))
## (kn.bic.air.m<-cub_spline_kn(air$xtab, air$ytab,
##   method="m", type="BIC"))
## (kn.bic.green.m<-cub_spline_kn(log(green$COST), 
##   log(green$OUTPUT), method="m", type="BIC"))
## (kn.bic.air.mc<-cub_spline_kn(air$xtab, air$ytab,
##   method="mc", type="BIC"))
## (kn.bic.green.mc<-cub_spline_kn(log(green$COST), 
##   log(green$OUTPUT), method="mc", type="BIC"))


###################################################
### code chunk number 30: ex-npbr.rnw:867-879 (eval = FALSE)
###################################################
## y.cub.air.u<-cub_spline_est(air$xtab, air$ytab, 
##  x.air, kn=kn.bic.air.u, method="u")
## y.cub.green.u<-cub_spline_est(log(green$COST), 
## log(green$OUTPUT),x.green,kn=kn.bic.green.u,method="u")
## y.cub.air.m<-cub_spline_est(air$xtab, air$ytab, 
##  x.air, kn=kn.bic.air.m, method="m")
## y.cub.green.m<-cub_spline_est(log(green$COST), 
## log(green$OUTPUT),x.green,kn=kn.bic.green.m,method="m") 
## y.cub.air.mc<-cub_spline_est(air$xtab, air$ytab, 
##  x.air, kn=kn.bic.air.mc, method="mc")
## y.cub.green.mc<-cub_spline_est(log(green$COST), 
## log(green$OUTPUT),x.green,kn=kn.bic.green.mc,method="mc") 


###################################################
### code chunk number 31: ex-npbr.rnw:883-899 (eval = FALSE)
###################################################
## plot(x.air, y.cub.air.u, lty=1, lwd=4, col="green", 
##  type="l", xlab="input", ylab="output")   
## lines(x.air, y.cub.air.m, lty=2, lwd=4, col="cyan")   
## lines(x.air, y.cub.air.mc, lty=3, lwd=4, col="magenta")   
## points(ytab~xtab, data=air)
## legend("topleft", col=c("green", "cyan","magenta"), 
##  lty=c(1,2,3), legend=c("unconstrained", "monotone",
##  "monotone+concave"), lwd=4, cex=0.8)
## plot(x.green, y.cub.green.u, lty=1, lwd=4, col="green", 
##  type="l", xlab="log(COST)", ylab="log(OUTPUT)")   
## lines(x.green, y.cub.green.m, lty=2, lwd=4,  col="cyan")   
## lines(x.green, y.cub.green.mc, lty=3, lwd=4,  col="magenta")   
## points(log(OUTPUT)~log(COST), data=green)
## legend("topleft", col=c("green","cyan","magenta"), 
## lty=c(1,2,3), legend=c("unconstrained", "monotone",  
##  "monotone+concave"), lwd=4, cex=0.8)


###################################################
### code chunk number 32: ex-npbr.rnw:903-923 (eval = FALSE)
###################################################
## .PngNo <- .PngNo + 1; name.file <- paste("Fig-bitmap-", .PngNo, ".pdf", sep="")
## pdf(file=name.file, width = 14, height = 7, pointsize = 14, bg = "white")
## op=par(mfrow=c(1,2),mar=c(3,3.1,2.1,2.1),mgp=c(2,.4,0),oma=c(0,0,0,0),cex.lab=1.2)
## plot(x.air, y.cub.air.u, lty=1, lwd=4, col="green", 
##  type="l", xlab="input", ylab="output")   
## lines(x.air, y.cub.air.m, lty=2, lwd=4, col="cyan")   
## lines(x.air, y.cub.air.mc, lty=3, lwd=4, col="magenta")   
## points(ytab~xtab, data=air)
## legend("topleft", col=c("green", "cyan","magenta"), 
##  lty=c(1,2,3), legend=c("unconstrained", "monotone",
##  "monotone+concave"), lwd=4, cex=0.8)
## plot(x.green, y.cub.green.u, lty=1, lwd=4, col="green", 
##  type="l", xlab="log(COST)", ylab="log(OUTPUT)")   
## lines(x.green, y.cub.green.m, lty=2, lwd=4,  col="cyan")   
## lines(x.green, y.cub.green.mc, lty=3, lwd=4,  col="magenta")   
## points(log(OUTPUT)~log(COST), data=green)
## legend("topleft", col=c("green","cyan","magenta"), 
## lty=c(1,2,3), legend=c("unconstrained", "monotone",  
##  "monotone+concave"), lwd=4, cex=0.8)
## par(op)
## dev.null <- dev.off()
## cat("\\includegraphics[width=0.9\\textwidth]{", name.file, "}\n\n", sep="")


###################################################
### code chunk number 33: ex-npbr.rnw:970-972 (eval = FALSE)
###################################################
## h.records.u<- loc_est_bw(records$year, 1/records$result, 
##  x.records, h=2, B=100, method="u")


###################################################
### code chunk number 34: ex-npbr.rnw:974-975 (eval = FALSE)
###################################################
## (h.records.u<-22.5)


###################################################
### code chunk number 35: ex-npbr.rnw:977-979 (eval = FALSE)
###################################################
## h.air.u<- loc_est_bw(air$xtab, air$ytab, x.air, 
##  h=2, B=100, method="u")


###################################################
### code chunk number 36: ex-npbr.rnw:981-982 (eval = FALSE)
###################################################
## (h.air.u<-3.612396)  


###################################################
### code chunk number 37: ex-npbr.rnw:984-986 (eval = FALSE)
###################################################
## h.air.m<- loc_est_bw(air$xtab, air$ytab, x.air,  
##  h=2, B=100, method="m") 


###################################################
### code chunk number 38: ex-npbr.rnw:988-989 (eval = FALSE)
###################################################
## (h.air.m<-3.638097)


###################################################
### code chunk number 39: ex-npbr.rnw:991-993 (eval = FALSE)
###################################################
## h.nucl.u <- loc_est_bw(nuclear$xtab, nuclear$ytab, 
##  x.nucl, h=40, B=100, method="u")


###################################################
### code chunk number 40: ex-npbr.rnw:995-996 (eval = FALSE)
###################################################
## (h.nucl.u<-79.11877)


###################################################
### code chunk number 41: ex-npbr.rnw:998-1000 (eval = FALSE)
###################################################
## h.nucl.m <- loc_est_bw(nuclear$xtab, nuclear$ytab, 
##  x.nucl, h=40, B=100, method="m")


###################################################
### code chunk number 42: ex-npbr.rnw:1002-1003 (eval = FALSE)
###################################################
## (h.nucl.m<-79.12)


###################################################
### code chunk number 43: ex-npbr.rnw:1007-1017 (eval = FALSE)
###################################################
## y.records.u<-loc_est(records$year, 1/records$result, 
##  x.records, h=h.records.u, method="u")
## y.air.u<-loc_est(air$xtab, air$ytab, x.air, h=h.air.u,
##  method="u")
## y.air.m<-loc_est(air$xtab, air$ytab, x.air, h=h.air.m, 
##  method="m")
## y.nucl.u<-loc_est(nuclear$xtab, nuclear$ytab, x.nucl, 
##  h=h.nucl.u, method="u")
## y.nucl.m<-loc_est(nuclear$xtab, nuclear$ytab, x.nucl, 
##  h=h.nucl.m, method="m") 


###################################################
### code chunk number 44: ex-npbr.rnw:1021-1038 (eval = FALSE)
###################################################
## plot(x.records, 1/y.records.u, lty=1, lwd=4, 
##  col="magenta", type="l")
## points(result~year, data=records)
## legend("topright",legend="unconstrained", col="magenta", 
##  lwd=4, lty=1)
##  
## plot(x.air, y.air.u, lty=1, lwd=4, col="magenta", type="l")
## lines(x.air, y.air.m, lty=2, lwd=4, col="cyan") 
## points(ytab~xtab, data=air)
## legend("topleft",legend=c("unconstrained", "improved"), 
##  col=c("magenta","cyan"), lwd=4, lty=c(1,2))
## 
## plot(x.nucl, y.nucl.u, lty=1, lwd=4, col="magenta", type="l")
## lines(x.nucl, y.nucl.m, lty=2, lwd=4, col="cyan") 
## points(ytab~xtab, data=nuclear)
## legend("topleft",legend=c("unconstrained", "improved"), 
##  col=c("magenta","cyan"), lwd=4, lty=c(1,2))


###################################################
### code chunk number 45: ex-npbr.rnw:1042-1061 (eval = FALSE)
###################################################
## .PngNo <- .PngNo + 1; name.file <- paste("Fig-bitmap-", .PngNo, ".pdf", sep="")
## pdf(file=name.file, width = 18, height = 7, pointsize = 14, bg = "white")
## op=par(mfrow=c(1,3),mar=c(3,3.1,2.1,2.1),mgp=c(2,.4,0),oma=c(0,0,0,0),cex.lab=1.2)
## plot(x.records, 1/y.records.u, lty=1, lwd=4, col="magenta", type="l",  xlab="year", ylab="1500m record")
## points(result~year, data=records)
## legend("topright",legend="unconstrained", col="magenta", lwd=4, lty=1)
##  
## plot(x.air, y.air.u, lty=1, lwd=4, col="magenta", type="l", xlab="input", ylab="output")
## lines(x.air, y.air.m, lty=2, lwd=4, col="cyan") 
## points(ytab~xtab, data=air)
## legend("topleft",legend=c("unconstrained", "improved"), 
##  col=c("magenta","cyan"), lwd=4, lty=c(1,2))
## 
## plot(x.nucl, y.nucl.u, lty=1, lwd=4, col="magenta", type="l", ylim=range(nuclear$ytab), xlab="temperature", ylab="toughness")
## lines(x.nucl, y.nucl.m, lty=2, lwd=4, col="cyan") 
## points(ytab~xtab, data=nuclear)
## legend("topleft",legend=c("unconstrained", "improved"), 
##  col=c("magenta","cyan"), lwd=4, lty=c(1,2))
##  par(op)
## dev.null <- dev.off()
## cat("\\includegraphics[width=0.9\\textwidth]{", name.file, "}\n\n", sep="")


###################################################
### code chunk number 46: ex-npbr.rnw:1099-1103 (eval = FALSE)
###################################################
## loc_max_1stage<-loc_max(log(green$COST), log(green$OUTPUT), 
##  x.green, h=0.5, type="one-stage")
## loc_max_2stage<-loc_max(log(green$COST), log(green$OUTPUT), 
##  x.green, h=0.5, type="two-stage")  


###################################################
### code chunk number 47: ex-npbr.rnw:1112-1116 (eval = FALSE)
###################################################
## require("np")
## bw <- npcdistbw(log(OUTPUT)~log(COST), data=green, 
##  cykertype = "uniform", bwtype="adaptive_nn")$xbw
## (h.opt<-max(bw, max(diff(sort(log(green$COST))))/2))


###################################################
### code chunk number 48: ex-npbr.rnw:1118-1119 (eval = FALSE)
###################################################
## (h.opt=0.4152283)


###################################################
### code chunk number 49: ex-npbr.rnw:1127-1131 (eval = FALSE)
###################################################
## loc_max_1stage.opt<-loc_max(log(green$COST), log(green$OUTPUT), 
##  x.green, h=h.opt, type="one-stage")
## loc_max_2stage.opt<-loc_max(log(green$COST), log(green$OUTPUT), 
##  x.green, h=h.opt, type="two-stage") 


###################################################
### code chunk number 50: ex-npbr.rnw:1134-1144 (eval = FALSE)
###################################################
## plot(log(OUTPUT)~log(COST), data=green)
## lines(x.green, loc_max_1stage, lty=1, col="magenta")
## lines(x.green, loc_max_2stage, lty=2, col="cyan")
## legend("topleft",legend=c("one-stage", "two-stage"), 
##  col=c("magenta","cyan"), lty=c(1,2))
## plot(log(OUTPUT)~log(COST), data=green)
## lines(x.green, loc_max_1stage.opt, lty=1, col="magenta")
## lines(x.green, loc_max_2stage.opt, lty=2, col="cyan")
## legend("topleft",legend=c("one-stage", "two-stage"), 
##  col=c("magenta","cyan"), lty=c(1,2)) 


###################################################
### code chunk number 51: ex-npbr.rnw:1149-1163 (eval = FALSE)
###################################################
## .PngNo <- .PngNo + 1; name.file <- paste("Fig-bitmap-", .PngNo, ".pdf", sep="")
## pdf(file=name.file, width = 14, height = 7, pointsize = 14, bg = "white")
## op=par(mfrow=c(1,2),mar=c(3,3.1,2.1,2.1),mgp=c(2,.4,0),oma=c(0,0,0,0),cex.lab=1.2)
## plot(log(OUTPUT)~log(COST), data=green, main="Peng and Gijbels choice")
## lines(x.green, loc_max_1stage, lty=1, lwd=2, col="magenta")
## lines(x.green, loc_max_2stage, lty=2, lwd=2, col="cyan")
## legend("topleft",legend=c("one-stage", "two-stage"), 
##  col=c("magenta","cyan"), lty=c(1,2),lwd=2)
## plot(log(OUTPUT)~log(COST), data=green, main="Automatic selection")
## lines(x.green, loc_max_1stage.opt, lty=1, lwd=2,col="magenta")
## lines(x.green, loc_max_2stage.opt, lty=2, lwd=2,col="cyan")
## legend("topleft",legend=c("one-stage", "two-stage"), 
##  col=c("magenta","cyan"), lty=c(1,2),lwd=2) 
## par(op)
## dev.null <- dev.off()
## cat("\\includegraphics[width=0.9\\textwidth]{", name.file, "}\n\n", sep="")


###################################################
### code chunk number 52: ex-npbr.rnw:1291-1294 (eval = FALSE)
###################################################
## require("npbr")
## require("np")
## data("green")


###################################################
### code chunk number 53: ex-npbr.rnw:1297-1302 (eval = FALSE)
###################################################
## require("np")
## data("green")
## (h.pr.green.m<-kern_smooth_bw(log(green$COST), 
##  log(green$OUTPUT), method="m", technique="pr",
##  bw_method="cv"))


###################################################
### code chunk number 54: ex-npbr.rnw:1304-1305 (eval = FALSE)
###################################################
## (h.pr.green.m<-0.8304566)


###################################################
### code chunk number 55: ex-npbr.rnw:1307-1310 (eval = FALSE)
###################################################
## (h.noh.green.m<-kern_smooth_bw(log(green$COST), 
##  log(green$OUTPUT), method="m", technique="noh",
##  bw_method="bic"))


###################################################
### code chunk number 56: ex-npbr.rnw:1315-1323 (eval = FALSE)
###################################################
## x.green <- seq(min(log(green$COST)), max(log(green$COST)),
##                length.out=101)
## y.pr.green.m<-kern_smooth(log(green$COST), 
##  log(green$OUTPUT), x.green, h=h.pr.green.m,
##  method="m", technique="pr")
## y.noh.green.m<-kern_smooth(log(green$COST), 
##  log(green$OUTPUT), x.green, h=h.noh.green.m,
##  method="m", technique="noh")


###################################################
### code chunk number 57: ex-npbr.rnw:1327-1333 (eval = FALSE)
###################################################
## plot(log(OUTPUT)~log(COST), data=green, xlab="log(COST)", 
##  ylab="log(OUTPUT)") 
## lines(x.green, y.pr.green.m, lwd=4, lty=3, col="red") 
## lines(x.green, y.noh.green.m, lwd=4, lty=3, col="blue")  
## legend("topleft", col=c("blue","red"), 
## lty=3, legend=c("noh","pr"), lwd=4, cex=0.8) 


###################################################
### code chunk number 58: ex-npbr.rnw:1337-1347 (eval = FALSE)
###################################################
## .PngNo <- .PngNo + 1; name.file <- paste("Fig-bitmap-", .PngNo, ".pdf", sep="")
## pdf(file=name.file, width = 14, height = 7, pointsize = 14, bg = "white")
## #op=par(mfrow=c(1,2),mar=c(3,3.1,2.1,2.1),mgp=c(2,.4,0),oma=c(0,0,0,0),cex.lab=1.2)
## plot(log(OUTPUT)~log(COST), data=green, xlab="log(COST)", 
##  ylab="log(OUTPUT)") 
## lines(x.green, y.pr.green.m, lwd=4, lty=3, col="red")  
## lines(x.green, y.noh.green.m, lwd=4, lty=3, col="blue") 
## legend("topleft", col=c("blue","red"), 
## lty=3, legend=c("noh","pr"), lwd=4, cex=0.8)
## #par(op)
## dev.null <- dev.off()
## cat("\\includegraphics[width=0.8\\textwidth]{", name.file, "}\n\n", sep="")


###################################################
### code chunk number 59: ex-npbr.rnw:1430-1432 (eval = FALSE)
###################################################
## x.post<- seq(post$xinput[100],max(post$xinput), 
##  length.out=100) 


###################################################
### code chunk number 60: ex-npbr.rnw:1435-1436 (eval = FALSE)
###################################################
## rho<-2


###################################################
### code chunk number 61: ex-npbr.rnw:1439-1441 (eval = FALSE)
###################################################
## best_kn.1<-kopt_momt_pick(post$xinput, post$yprod, 
##  x.post, rho=rho)


###################################################
### code chunk number 62: ex-npbr.rnw:1444-1446 (eval = FALSE)
###################################################
## rho_momt<-rho_momt_pick(post$xinput, post$yprod, 
##  x.post, method="moment")


###################################################
### code chunk number 63: ex-npbr.rnw:1449-1451 (eval = FALSE)
###################################################
## best_kn.2<-kopt_momt_pick(post$xinput, post$yprod,
##   x.post, rho=rho_momt) 


###################################################
### code chunk number 64: ex-npbr.rnw:1453-1469 (eval = FALSE)
###################################################
## rho_momt<-c(1.993711,2.360920,2.245450,3.770526,2.724960,3.667846,4.026203,
## 2.281109,1.363260,1.150343,2.567832,2.228400,3.106491,2.592477,
## 2.233479,2.040209,1.916878,1.494831,1.961430,1.930942,1.927990,
## 1.833530,1.808632,1.758135,1.717626,1.686540,1.707200,1.711357,
## 1.720839,1.704845,1.678985,1.686872,1.686907,1.747732,1.741290,
## 1.792388,1.805144,1.855829,1.919817,1.929348,2.046588,2.135351,
## 2.196834,2.224797,2.221043,2.290578,2.390179,2.042884,2.087287,
## 2.158198,2.173314,2.260872,2.311427,1.865147,1.874019,1.913673,
## 1.922869,1.918484,1.949220,1.961016,1.998101,2.023605,2.041663,
## 2.067775,2.088982,2.107949,2.152688,2.170959,1.283350,1.285458,
## 1.295437,1.296902,1.316896,1.331668,1.330163,1.339701,1.322501,
## 1.326488,1.373837,1.392537,1.419458,1.426513,1.448544,1.473716,
## 1.517720,1.549229,1.561259,1.567216,1.580512,1.647293,1.672556,
## 1.750994,1.743083,1.801643,1.823678,1.869798,1.906898,1.873269,
## 1.893699,1.916469)
## best_kn.2<-kopt_momt_pick(post$xinput, post$yprod, x.post, rho=rho_momt)


###################################################
### code chunk number 65: ex-npbr.rnw:1476-1479 (eval = FALSE)
###################################################
## rho_trimmean<-mean(rho_momt, trim=0.00)
## best_kn.3<-kopt_momt_pick(post$xinput, post$yprod,
##   x.post, rho=rho_trimmean) 


###################################################
### code chunk number 66: ex-npbr.rnw:1482-1488 (eval = FALSE)
###################################################
## res.momt.1<-dfs_momt(post$xinput, post$yprod, x.post, 
##  rho=rho, k=best_kn.1)
## res.momt.2<-dfs_momt(post$xinput, post$yprod, x.post, 
##  rho=rho_momt, k=best_kn.2)
## res.momt.3<-dfs_momt(post$xinput, post$yprod, x.post, 
##  rho=rho_trimmean, k=best_kn.3) 


###################################################
### code chunk number 67: ex-npbr.rnw:1491-1506 (eval = FALSE)
###################################################
## plot(yprod~xinput, data=post, xlab="Quantity of labor", 
##  ylab="Volume of delivered mail")
## lines(x.post, res.momt.1[,1], lty=1, col="cyan")  
## lines(x.post, res.momt.1[,2], lty=3, col="magenta")  
## lines(x.post, res.momt.1[,3], lty=3, col="magenta")  
## plot(yprod~xinput, data=post, xlab="Quantity of labor", 
##  ylab="Volume of delivered mail")
## lines(x.post, res.momt.2[,1], lty=1, col="cyan")  
## lines(x.post, res.momt.2[,2], lty=3, col="magenta")  
## lines(x.post, res.momt.2[,3], lty=3, col="magenta") 
## plot(yprod~xinput, data=post, xlab="Quantity of labor", 
##  ylab="Volume of delivered mail")
## lines(x.post, res.momt.3[,1], lty=1, col="cyan")  
## lines(x.post, res.momt.3[,2], lty=3, col="magenta")  
## lines(x.post, res.momt.3[,3], lty=3, col="magenta")  


###################################################
### code chunk number 68: ex-npbr.rnw:1510-1526 (eval = FALSE)
###################################################
## .PngNo <- .PngNo + 1; name.file <- paste("Fig-bitmap-", .PngNo, ".pdf", sep="")
## pdf(file=name.file, width = 18, height = 7, pointsize = 14, bg = "white")
## op=par(mfrow=c(1,3),mar=c(3,3.1,2.1,2.1),mgp=c(2,.4,0),oma=c(0,0,0,0),cex.lab=1.2)
## plot(yprod~xinput, data=post,  col="grey", xlab="Quantity of labor", ylab="Volume of delivered mail")
## lines(x.post, res.momt.1[,1], lty=1, lwd=2, col="cyan")  
## lines(x.post, res.momt.1[,2], lty=3, lwd=4, col="magenta")  
## lines(x.post, res.momt.1[,3], lty=3, lwd=4, col="magenta")  
## plot(yprod~xinput, data=post,  col="grey", xlab="Quantity of labor", ylab="Volume of delivered mail")
## lines(x.post, res.momt.2[,1], lty=1, lwd=2, col="cyan")  
## lines(x.post, res.momt.2[,2], lty=3, lwd=4, col="magenta")  
## lines(x.post, res.momt.2[,3], lty=3, lwd=4, col="magenta")
## plot(yprod~xinput, data=post,  col="grey", xlab="Quantity of labor", ylab="Volume of delivered mail")
## lines(x.post, res.momt.3[,1], lty=1, lwd=2, col="cyan")  
## lines(x.post, res.momt.3[,2], lty=3, lwd=4, col="magenta")  
## lines(x.post, res.momt.3[,3], lty=3, lwd=4, col="magenta")    
## par(op)
## dev.null <- dev.off()
## cat("\\includegraphics[width=0.9\\textwidth]{", name.file, "}\n\n", sep="")


###################################################
### code chunk number 69: ex-npbr.rnw:1581-1582 (eval = FALSE)
###################################################
## rho<-2


###################################################
### code chunk number 70: ex-npbr.rnw:1585-1587 (eval = FALSE)
###################################################
## best_kn.1<-kopt_momt_pick(post$xinput, post$yprod, 
##  x.post, method="pickands", rho=rho)


###################################################
### code chunk number 71: ex-npbr.rnw:1590-1592 (eval = FALSE)
###################################################
## rho_pick<-rho_momt_pick(post$xinput, post$yprod, 
##  x.post, method="pickands")


###################################################
### code chunk number 72: ex-npbr.rnw:1595-1597 (eval = FALSE)
###################################################
## best_kn.2<-kopt_momt_pick(post$xinput, post$yprod,
##   x.post, method="pickands", rho=rho_pick)


###################################################
### code chunk number 73: ex-npbr.rnw:1601-1604 (eval = FALSE)
###################################################
## rho_trimmean<-mean(rho_pick, trim=0.00)
## best_kn.3<-kopt_momt_pick(post$xinput, post$yprod,
##   x.post, rho=rho_trimmean, method="pickands") 


###################################################
### code chunk number 74: ex-npbr.rnw:1607-1613 (eval = FALSE)
###################################################
## res.pick.1<-dfs_pick(post$xinput, post$yprod, x.post, 
##  rho=rho, k=best_kn.1)
## res.pick.2<-dfs_pick(post$xinput, post$yprod, x.post, 
##  rho=rho_pick, k=best_kn.2)
## res.pick.3<-dfs_pick(post$xinput, post$yprod, x.post, 
##  rho=rho_trimmean, k=best_kn.3) 


###################################################
### code chunk number 75: ex-npbr.rnw:1617-1632 (eval = FALSE)
###################################################
## plot(yprod~xinput, data=post, xlab="Quantity of labor", 
##  ylab="Volume of delivered mail")
## lines(x.post, res.pick.1[,1], lty=1, col="cyan")  
## lines(x.post, res.pick.1[,2], lty=3, col="magenta")  
## lines(x.post, res.pick.1[,3], lty=3, col="magenta")  
## plot(yprod~xinput, data=post, xlab="Quantity of labor", 
##  ylab="Volume of delivered mail")
## lines(x.post, res.pick.2[,1], lty=1, col="cyan")  
## lines(x.post, res.pick.2[,2], lty=3, col="magenta")  
## lines(x.post, res.pick.2[,3], lty=3, col="magenta")  
## plot(yprod~xinput, data=post, xlab="Quantity of labor", 
##  ylab="Volume of delivered mail")
## lines(x.post, res.pick.3[,1], lty=1, col="cyan")  
## lines(x.post, res.pick.3[,2], lty=3, col="magenta")  
## lines(x.post, res.pick.3[,3], lty=3, col="magenta")  


###################################################
### code chunk number 76: ex-npbr.rnw:1636-1652 (eval = FALSE)
###################################################
## .PngNo <- .PngNo + 1; name.file <- paste("Fig-bitmap-", .PngNo, ".pdf", sep="")
## pdf(file=name.file, width = 18, height = 7, pointsize = 14, bg = "white")
## op=par(mar=c(3,3.1,2.1,2.1),mgp=c(2,.4,0),oma=c(0,0,0,0),cex.lab=1.2, mfrow=c(1,3))
## plot(yprod~xinput, data=post, col="grey", xlab="Quantity of labor", ylab="Volume of delivered mail")
## lines(x.post, res.pick.1[,1], lty=1, lwd=2, col="cyan")  
## lines(x.post, res.pick.1[,2], lty=3, lwd=4, col="magenta")  
## lines(x.post, res.pick.1[,3], lty=3, lwd=4, col="magenta")  
## plot(yprod~xinput, data=post, col="grey", xlab="Quantity of labor", ylab="Volume of delivered mail")
## lines(x.post, res.pick.2[,1], lty=1, lwd=2, col="cyan")  
## lines(x.post, res.pick.2[,2], lty=3, lwd=4, col="magenta")  
## lines(x.post, res.pick.2[,3], lty=3, lwd=4, col="magenta")  
## plot(yprod~xinput, data=post, col="grey", xlab="Quantity of labor", ylab="Volume of delivered mail")
## lines(x.post, res.pick.3[,1], lty=1, lwd=2, col="cyan")  
## lines(x.post, res.pick.3[,2], lty=3, lwd=4, col="magenta")  
## lines(x.post, res.pick.3[,3], lty=3, lwd=4, col="magenta") 
## par(op)
## dev.null <- dev.off()
## cat("\\includegraphics[width=0.9\\textwidth]{", name.file, "}\n\n", sep="")


###################################################
### code chunk number 77: ex-npbr.rnw:1772-1773 (eval = FALSE)
###################################################
## rho<-2


###################################################
### code chunk number 78: ex-npbr.rnw:1776-1780 (eval = FALSE)
###################################################
## best_cm.1<- mopt_pwm(post$xinput, post$yprod, 
##  x.post, a=2, rho=rho, wind.coef=0.1)
## res.pwm.1<- dfs_pwm(post$xinput, post$yprod, x.post, 
##  coefm=best_cm.1, a=2, rho=rho) 


###################################################
### code chunk number 79: ex-npbr.rnw:1782-1807 (eval = FALSE)
###################################################
## res.pwm.1<-matrix(c(3693.689,3391.852,3995.527,3698.156,3393.867,4002.446,3664.170,3378.494,3949.845,4052.185,3678.129,
## 4426.240,4583.859,4082.118,5085.599,4544.788,4058.775,5030.800,4514.398,4040.244,4988.552,4370.231,3966.160,4774.301,
## 4299.334,3920.889,4677.779,5996.560,5104.971,6888.149,5767.108,5006.010,6528.207,8134.150,6676.470,9591.829,7714.785,
## 6449.040,8980.529,7293.492,6169.318,8417.665,7182.106,6156.329,8207.884,6996.592,6065.089,7928.095,7237.842,6381.361,
## 8094.322,7028.955,6258.911,7798.999,6880.657,6164.981,7596.333,6871.169,6193.560,7548.777,7087.487,6442.267,7732.707,
## 7160.956,6542.813,7779.100,7473.784,6873.938,8073.630,7468.800,6889.968,8047.631,7423.528,6855.837,7991.220,7399.736,
## 6846.585,7952.887,7466.078,6922.605,8009.552,7695.183,7152.966,8237.401,7674.784,7146.509,8203.058,7632.812,7115.122,
## 8150.503,7656.791,7143.767,8169.815,7694.408,7188.657,8200.158,7743.059,7245.345,8240.773,7784.993,7295.642,8274.344,
## 7768.774,7282.830,8254.718,7744.851,7268.383,8221.318,7737.486,7263.321,8211.652,7816.011,7351.348,8280.673,8050.885,
## 7582.127,8519.643,8065.387,7600.336,8530.437,8379.410,7909.951,8848.869,8691.761,8217.325,9166.196,8697.871,8234.138,
## 9161.604,8671.017,8213.569,9128.464,8673.969,8215.666,9132.272,9115.328,8640.232,9590.424,9388.969,8909.159,9868.779,
## 9496.951,9037.316,9956.587,9500.138,9050.094,9950.183,9709.822,9256.602,10163.042,9684.555,9236.256,10132.853,9973.046,
## 9498.580,10447.511,10399.506,9868.298,10930.713,10379.734,9851.823,10907.644,10406.428,9882.258,10930.599,10473.119,
## 9956.109,10990.130,10524.562,10009.960,11039.163,10520.942,10006.962,11034.921,10583.200,10072.521,11093.879,10643.052,
## 10134.391,11151.712,10638.578,10141.266,11135.889,10625.220,10133.219,11117.221,10611.912,10122.349,11101.475,10731.387,
## 10240.820,11221.954,10709.384,10223.295,11195.474,10697.309,10217.597,11177.020,11335.434,10745.236,11925.632,11320.899,
## 10733.487,11908.310,11357.880,10778.002,11937.758,11364.218,10782.596,11945.840,11345.092,10767.871,11922.313,11340.068,
## 10763.661,11916.475,11483.711,10903.637,12063.785,11543.050,10967.724,12118.376,11551.772,10983.171,12120.374,11576.991,
## 11015.011,12138.971,11546.826,10990.595,12103.056,11536.971,10982.598,12091.344,11880.266,11316.796,12443.736,11950.866,
## 11394.222,12507.509,11994.852,11452.743,12536.962,11989.201,11447.978,12530.425,11964.985,11431.400,12498.571,12014.366,
## 11490.019,12538.714,12291.219,11762.119,12820.319,12390.599,11866.754,12914.443,12356.355,11845.739,12866.971,12336.726,
## 11830.709,12842.744,12581.846,12055.270,13108.422,12996.801,12447.918,13545.684,12976.842,12437.372,13516.311,13558.110,
## 12982.361,14133.860,13541.724,12968.784,14114.664,13505.548,12946.694,14064.402,13490.559,12933.933,14047.185,13622.065,
## 13062.711,14181.420,13611.193,13061.849,14160.538,13586.358,13044.210,14128.507,13570.022,13030.314,14109.731,13564.974,
## 13029.502,14100.446),100,3,byrow=TRUE)


###################################################
### code chunk number 80: ex-npbr.rnw:1810-1813 (eval = FALSE)
###################################################
## rho_pwm<-rho_pwm(post$xinput, post$yprod, 
##  x.post, a=2, lrho=1, urho=Inf)
## rho_pwm_trim<-mean(rho_pwm, trim=0.00)


###################################################
### code chunk number 81: ex-npbr.rnw:1815-1829 (eval = FALSE)
###################################################
## rho_pwm<-c(1.023594,1.024185,1.039690,1.052159,1.024773,1.039298,1.927103,
## 1.837867,1.677647,1.550235,1.454431,1.379524,1.310404,1.260026,1.234148,
## 1.203759,1.178933,1.170742,1.161470,1.149357,1.198126,1.314633,1.515514,
## 1.599483,1.665202,1.722867,1.817279,1.879157,1.945815,1.990754,2.031745,
## 2.105846,2.165757,2.237778,2.259292,2.295147,2.315305,2.397803,2.079720,
## 2.092493,2.206906,2.359096,2.385121,2.400971,2.422183,2.024605,2.243093,
## 2.294995,2.310249,2.317997,2.334523,1.181741,1.189079,1.191125,1.199886,
## 1.204921,1.205482,1.208576,1.202842,1.209828,1.215217,1.195819,1.198989,
## 1.185889,1.188102,1.190547,1.580017,1.581408,1.586728,1.589427,1.590941,
## 1.592341,1.706439,1.721817,1.724899,1.729585,1.732874,1.735634,2.072950,
## 2.097440,2.137919,2.140791,2.144488,2.162599,5.093260,6.434999,6.410624,
## 6.412833,4.358651,3.123234,3.148811,1.078259,1.079294,1.074666,1.076102,
## 1.082311,1.083827,1.077977,1.079619,1.080355)
## rho_pwm_trim<-mean(rho_pwm, trim=0.00)


###################################################
### code chunk number 82: ex-npbr.rnw:1832-1840 (eval = FALSE)
###################################################
## best_cm.2<- mopt_pwm(post$xinput, post$yprod, 
##  x.post, a=2, rho = rho_pwm)
## best_cm.3<- mopt_pwm(post$xinput, post$yprod, 
##  x.post, a=2, rho = rho_pwm_trim)
## res.pwm.2<- dfs_pwm(post$xinput, post$yprod, 
##  x.post, coefm=best_cm.2, rho=rho_pwm)
## res.pwm.3<- dfs_pwm(post$xinput, post$yprod, 
##  x.post, coefm=best_cm.3, rho=rho_pwm_trim)


###################################################
### code chunk number 83: ex-npbr.rnw:1842-1910 (eval = FALSE)
###################################################
## res.pwm.2<-matrix(c(3423.634,3133.666,3713.601,3420.217,3127.762,3712.672,
## 3403.614,3117.938,3689.290,3709.516,3335.461,4083.571,4093.768,3592.028,4595.509,
## 4056.693,3570.681,4542.706,4476.577,4002.423,4950.731,4294.037,3889.966,4698.107,
## 4152.061,3773.616,4530.506,5570.299,4678.710,6461.887,5307.500,4546.401,6068.599,
## 7158.666,5700.986,8616.345,6731.634,5465.889,7997.379,6339.216,5215.043,7463.390,
## 6232.965,5207.188,7258.743,6068.923,5137.420,7000.426,6265.447,5408.967,7121.928,
## 6137.926,5367.883,6907.970,6035.611,5319.935,6751.287,6039.558,5361.950,6717.166,
## 6290.624,5645.404,6935.844,6477.025,5858.882,7095.168,6966.139,6366.293,7565.985,
## 7056.248,6477.417,7635.079,7083.572,6515.880,7651.263,7125.345,6572.194,7678.496,
## 7285.332,6741.859,7828.805,7571.149,7028.932,8113.367,7620.769,7092.495,8149.044,
## 7623.786,7106.095,8141.477,7687.649,7174.624,8200.673,7797.465,7291.714,8303.215,
## 7904.813,7407.099,8402.527,8013.950,7524.599,8503.301,8016.589,7530.645,8502.533,
## 8020.621,7544.154,8497.088,8028.042,7553.876,8502.208,8182.375,7717.713,8647.037,
## 8127.045,7658.287,8595.803,8153.546,7688.495,8618.597,8586.068,8116.609,9055.527,
## 9069.230,8594.794,9543.665,9094.677,8630.944,9558.409,9076.375,8618.927,9533.822,
## 9100.556,8642.253,9558.859,9142.547,8667.451,9617.643,9669.278,9189.468,10149.088,
## 9830.317,9370.682,10289.953,9842.657,9392.612,10292.701,10069.847,9616.627,10523.067,
## 10058.510,9610.212,10506.808,9019.517,8545.051,9493.982,9359.357,8828.149,9890.564,
## 9351.436,8823.526,9879.347,9389.890,8865.720,9914.060,9469.730,8952.719,9986.740,
## 9517.669,9003.067,10032.270,9518.444,9004.464,10032.424,9571.469,9060.790,10082.148,
## 9634.284,9125.624,10142.945,9656.974,9159.662,10154.285,9631.020,9139.019,10123.022,
## 9627.760,9138.197,10117.324,9716.570,9226.003,10207.137,9709.098,9223.009,10195.188,
## 9713.317,9233.605,10193.029,10746.776,10156.577,11336.974,10738.178,10150.766,11325.589,
## 10786.367,10206.489,11366.245,10797.398,10215.776,11379.019,10786.170,10208.949,11363.391,
## 10783.844,10207.437,11360.252,11073.602,10493.528,11653.676,11154.008,10578.682,11729.334,
## 11170.647,10602.046,11739.248,11205.198,10643.218,11767.177,11184.990,10628.759,11741.220,
## 11180.277,10625.904,11734.650,11982.912,11419.442,12546.382,12087.762,11531.119,12644.406,
## 12185.205,11643.096,12727.315,12183.219,11641.995,12724.442,12159.491,11625.905,12693.077,
## 12231.392,11707.045,12755.740,16571.616,16042.516,17100.716,18554.292,18030.448,19078.137,
## 18317.125,17806.509,18827.741,18235.602,17729.584,18741.619,15878.405,15351.829,16404.981,
## 14649.501,14100.618,15198.383,14637.705,14098.236,15177.175,12136.540,11560.790,12712.289,
## 12130.230,11557.290,12703.170,12131.561,11572.707,12690.416,12124.857,11568.231,12681.483,
## 12246.209,11686.855,12805.564,12260.300,11710.956,12809.644,12246.975,11704.827,12789.123,
## 12240.520,11700.812,12780.229,12245.702,11710.229,12781.174),100,3,byrow=TRUE)
## res.pwm.3<-matrix(c(3646.226,3344.389,3948.064,3648.807,3344.517,3953.097,
## 3622.910,3337.234,3908.586,3997.221,3623.165,4371.276,4507.417,4005.677,5009.158,
## 4467.527,3981.515,4953.540,4436.192,3962.038,4910.346,4299.358,3895.287,4703.428,
## 4230.362,3851.917,4608.807,5853.342,4961.754,6744.931,5639.696,4878.597,6400.795,
## 7896.196,6438.517,9353.876,7498.835,6233.090,8764.580,7098.036,5973.863,8222.210,
## 6994.209,5968.432,8019.987,6819.883,5888.380,7751.386,7058.151,6201.671,7914.632,
## 6865.907,6095.864,7635.951,6727.713,6012.037,7443.389,6722.773,6045.165,7400.381,
## 6936.750,6291.530,7581.970,7009.811,6391.668,7627.954,7315.399,6715.554,7915.245,
## 7313.206,6734.375,7892.037,7270.222,6702.530,7837.913,7250.309,6697.159,7803.460,
## 7316.880,6773.406,7860.353,7540.428,6998.210,8082.645,7524.539,6996.264,8052.813,
## 7485.710,6968.019,8003.401,7510.348,6997.323,8023.372,7547.774,7042.024,8053.525,
## 7596.133,7098.420,8093.847,7640.059,7150.708,8129.410,7624.931,7138.987,8110.875,
## 7604.246,7127.779,8080.713,7598.825,7124.659,8072.991,7677.471,7212.809,8142.133,
## 7906.991,7438.233,8375.749,7921.832,7456.781,8386.883,8229.054,7759.595,8698.513,
## 8533.614,8059.179,9008.049,8542.872,8079.139,9006.605,8518.945,8061.497,8976.392,
## 8521.985,8063.682,8980.288,8948.660,8473.564,9423.756,9215.412,8735.602,9695.222,
## 9326.893,8867.258,9786.529,9334.011,8883.967,9784.056,9539.464,9086.243,9992.684,
## 9516.357,9068.059,9964.655,9796.244,9321.778,10270.709,10204.918,9673.710,10736.125,
## 10186.882,9658.972,10714.792,10213.718,9689.548,10737.888,10281.710,9764.700,10798.721,
## 10332.350,9817.748,10846.952,10328.829,9814.849,10842.809,10390.698,9880.019,10901.376,
## 10449.435,9940.774,10958.095,10448.894,9951.583,10946.206,10437.687,9945.686,10929.688,
## 10425.548,9935.984,10915.111,10542.273,10051.706,11032.840,10522.477,10036.387,11008.566,
## 10512.897,10033.185,10992.608,11123.680,10533.481,11713.878,11110.586,10523.175,11697.997,
## 11148.966,10569.088,11728.844,11155.662,10574.040,11737.283,11138.683,10561.462,11715.904,
## 11133.952,10557.545,11710.360,11272.850,10692.776,11852.924,11331.986,10756.659,11907.312,
## 11342.691,10774.090,11911.292,11369.501,10807.521,11931.481,11342.411,10786.180,11898.641,
## 11333.361,10778.988,11887.734,11668.327,11104.857,12231.797,11739.272,11182.628,12295.915,
## 11787.024,11244.914,12329.133,11781.696,11240.473,12322.920,11762.284,11228.699,12295.870,
## 11813.403,11289.055,12337.750,12083.976,11554.877,12613.076,12182.628,11658.783,12706.472,
## 12154.117,11643.501,12664.733,12136.689,11630.671,12642.707,12372.389,11845.813,12898.965,
## 12775.894,12227.011,13324.777,12759.797,12220.327,13299.266,13323.740,12747.991,13899.490,
## 13308.758,12735.818,13881.698,13279.889,12721.034,13838.743,13265.917,12709.291,13822.544,
## 13394.246,12834.892,13953.601,13387.143,12837.799,13936.488,13365.605,12823.457,13907.753,
## 13350.512,12810.804,13890.221,13346.982,12811.509,13882.454),100,3,byrow=TRUE)


###################################################
### code chunk number 84: ex-npbr.rnw:1913-1928 (eval = FALSE)
###################################################
## plot(yprod~xinput, data=post, xlab="Quantity of labor", 
##  ylab="Volume of delivered mail")
## lines(x.post, res.pwm.1[,1], lty=1, col="cyan")  
## lines(x.post, res.pwm.1[,2], lty=3, col="magenta")  
## lines(x.post, res.pwm.1[,3], lty=3, col="magenta")  
## plot(yprod~xinput, data=post, xlab="Quantity of labor", 
##  ylab="Volume of delivered mail")
## lines(x.post, res.pwm.2[,1], lty=1, col="cyan")  
## lines(x.post, res.pwm.2[,2], lty=3, col="magenta")  
## lines(x.post, res.pwm.2[,3], lty=3, col="magenta")
## plot(yprod~xinput, data=post, xlab="Quantity of labor", 
##  ylab="Volume of delivered mail")
## lines(x.post, res.pwm.3[,1], lty=1, col="cyan")  
## lines(x.post, res.pwm.3[,2], lty=3, col="magenta")  
## lines(x.post, res.pwm.3[,3], lty=3, col="magenta")   


###################################################
### code chunk number 85: ex-npbr.rnw:1932-1949 (eval = FALSE)
###################################################
## .PngNo <- .PngNo + 1; name.file <- paste("Fig-bitmap-", .PngNo, ".pdf", sep="")
## pdf(file=name.file, width = 18, height = 7, pointsize = 14, bg = "white")
## op=par(mar=c(3,3.1,2.1,2.1),mgp=c(2,.4,0),oma=c(0,0,0,0),cex.lab=1.2, mfrow=c(1,3))
## plot(yprod~xinput, data=post, col="grey", xlab="Quantity of labor", ylab="Volume of delivered mail")
## lines(x.post, res.pwm.1[,1], lty=1, lwd=2, col="cyan")  
## lines(x.post, res.pwm.1[,2], lty=3, lwd=4, col="magenta")  
## lines(x.post, res.pwm.1[,3], lty=3, lwd=4, col="magenta")  
## plot(yprod~xinput, data=post, col="grey", xlab="Quantity of labor", ylab="Volume of delivered mail")
## lines(x.post, res.pwm.2[,1], lty=1, lwd=2, col="cyan")  
## lines(x.post, res.pwm.2[,2], lty=3, lwd=4, col="magenta")  
## lines(x.post, res.pwm.2[,3], lty=3, lwd=4, col="magenta")  
## plot(yprod~xinput, data=post, xlab="Quantity of labor", col="grey", 
##  ylab="Volume of delivered mail")
## lines(x.post, res.pwm.3[,1], lty=1, lwd=2, col="cyan")  
## lines(x.post, res.pwm.3[,2], lty=3, lwd=4, col="magenta")  
## lines(x.post, res.pwm.3[,3], lty=3, lwd=4, col="magenta")  
## par(op)
## dev.null <- dev.off()
## cat("\\includegraphics[width=0.9\\textwidth]{", name.file, "}\n\n", sep="")


###################################################
### code chunk number 86: ex-npbr.rnw:2013-2014 (eval = FALSE)
###################################################
## require("npbr")


###################################################
### code chunk number 87: ex-npbr.rnw:2017-2038 (eval = FALSE)
###################################################
## N<-5
## x.sim <- seq(0, 1, length.out=1000)
## y.dea<-matrix(0, N, 1000)
## y.cub<-matrix(0, N, 1000)
## 
## Fron<-function(x) sqrt(x)
## 
## for (k in 1:N)
## { 
##   n=100; betav=0.5
##   xtab <- runif(n, 0, 1)
##   V <-rbeta(n, betav, betav)
##   ytab <- Fron(xtab)*V
##   cind<-which((x.sim>=min(xtab))&(x.sim<=max(xtab)))
##   x<-x.sim[cind]  
##   y.dea[k,cind]<-dea_est(xtab, ytab, x, type="dea")
##   kopt<-cub_spline_kn(xtab,ytab,method="mc",krange=1:20,
##                       type="BIC")
##   y.cub[k,cind]<-cub_spline_est(xtab,ytab,x,kn=kopt,
##                       method="mc",all.dea=FALSE)
## }


###################################################
### code chunk number 88: ex-npbr.rnw:2041-2065 (eval = FALSE)
###################################################
## require("npbr")
## evaluation<-function(MAT,xeval,true_vec)
## {
##   # internal function
##   denzero<-function(vec)
##   {
##     return(sum(vec!=0))
##   }
##   
##   nzr<-apply(MAT,1,denzero) 
##   nzc<-apply(MAT,2,denzero)
##   nzc_ind<-which(apply(MAT,2,denzero)!=0)
##   nz_mat<-matrix(as.numeric(MAT!=0),dim(MAT)[1],length(xeval),byrow=FALSE)
##   cmean<-rep(0,dim(MAT)[2])
##   temp<-apply(MAT,2,sum)
##   cmean[nzc_ind]<-temp[nzc_ind]*(1/nzc[nzc_ind])
##   
##   temp2<-apply((MAT-rep(1,dim(MAT)[1]) %*% t(cmean))^2 * nz_mat,2,sum)
##   IVAR<-mean(temp2[nzc_ind]*(1/nzc[nzc_ind]))
##   temp3<-(true_vec-cmean)^2
##   IBIAS<-mean(temp3[nzc_ind])
##   IMSE<-IBIAS+IVAR
##   return(list(IBIAS2=IBIAS,IVAR=IVAR,MISE=IMSE))
## }


###################################################
### code chunk number 89: ex-npbr.rnw:2068-2071 (eval = FALSE)
###################################################
## result.dea<-evaluation(y.dea,x.sim,Fron(x.sim))
## result.cub<-evaluation(y.cub,x.sim,Fron(x.sim))
## (cbind(result.dea,result.cub))


