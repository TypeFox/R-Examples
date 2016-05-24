### R code from vignette source 'nsRFA_ex02.Rnw'

###################################################
### code chunk number 1: nsRFA_ex02.Rnw:38-39
###################################################
library(nsRFA)


###################################################
### code chunk number 2: nsRFA_ex02.Rnw:42-43
###################################################
data(FEH1000)


###################################################
### code chunk number 3: nsRFA_ex02.Rnw:46-48 (eval = FALSE)
###################################################
## ls()
## help(FEH1000)


###################################################
### code chunk number 4: nsRFA_ex02.Rnw:56-65
###################################################
n <- tapply(am[,4],am[,1],length)
urbext <- cd[,"urbext1990"]
area <- cd[,"dtm_area"]
cd696 <- cd[(!is.nan(cd[,"dtm_area"]))&(!is.nan(cd[,"saar"]))&
            (!is.nan(cd[,"bfihost"]))&(n>7)&(urbext<0.025)&(area>0.5),]

fac <- factor(am[,"number"],levels=cd696[,"number"])
am696 <- am[!is.na(fac),]
#nlevels(as.factor(am696[,"number"]))


###################################################
### code chunk number 5: nsRFA_ex02.Rnw:74-78
###################################################
layout(matrix(c(1,2,3), 1, 3))
plot(cd696[c("dtm_area","saar")], pch=".", cex=2, log="x")
plot(cd696[c("dtm_area","bfihost")], pch=".", cex=2, log="x")
plot(cd696[c("saar","bfihost")], pch=".", cex=2)


###################################################
### code chunk number 6: nsRFA_ex02.Rnw:91-93
###################################################
Lmomenti696 <- t(sapply(split(am696[,4],am696[,1]),Lmoments))
Di <- discordancy(am696[,"am"], am696[,"number"])


###################################################
### code chunk number 7: nsRFA_ex02.Rnw:99-105
###################################################
par(mfrow=c(1,2))
 plot(Lmomenti696[,c("lca","lcv")],xlab="L-CA",ylab="L-CV",pch=".",cex=2); grid()
 points(Lmomenti696[(Di>3),c("lca","lcv")],pch=19,cex=.7)
 plot(Lmomenti696[,c("lca","lkur")],xlab="L-CA",ylab="L-kur",pch=".",cex=2); grid()
 points(Lmomenti696[(Di>3),c("lca","lkur")],pch=19,cex=.7)
par(mfrow=c(1,1))


###################################################
### code chunk number 8: nsRFA_ex02.Rnw:120-123
###################################################
sd(log(cd696[,"dtm_area"])) # 1.345515 (vs 1.34)
sd(log(cd696[,"saar"]))     # 0.38534 (vs 0.38)
sd(cd696[,"bfihost"])       # 0.1485239 (vs 0.15)


###################################################
### code chunk number 9: nsRFA_ex02.Rnw:125-133
###################################################
AREAterm <- log(cd696[,"dtm_area"])/(sd(log(cd696[,"dtm_area"]))*sqrt(2))
SAARterm <- log(cd696[,"saar"])/sd(log(cd696[,"saar"]))
BFIHOSTterm <- cd696[,"bfihost"]/sd(cd696[,"bfihost"])

distFEH <- dist(cbind(AREAterm,SAARterm,BFIHOSTterm))

roi.cd <- data.frame(cbind(AREAterm,SAARterm,BFIHOSTterm))
row.names(roi.cd) <- cd696[,"number"]


###################################################
### code chunk number 10: nsRFA_ex02.Rnw:136-144 (eval = FALSE)
###################################################
## roi01.50year <- new.env()
## for(i in 1:696) {
##  print(paste(i,"/ 696"))
##  assign(as.character(row.names(roi.cd)[i]), roi.st.year(roi.cd[i,],as.data.frame(roi.cd),
##      row.names(roi.cd),am696[,"am"],am696[,"number"],test="HW",station.year=250,Nsim=100), 
##      env=roi01.50year)
## }
## roi01.50year <- as.list(roi01.50year)


###################################################
### code chunk number 11: nsRFA_ex02.Rnw:146-148
###################################################
estrai.region <- function (x) {x$region}
estrai.test <- function (x) {x$test}


###################################################
### code chunk number 12: nsRFA_ex02.Rnw:150-156 (eval = FALSE)
###################################################
## regioni.50year <- sapply(roi01.50year, estrai.region)
## test.50year <- sapply(roi01.50year, estrai.test)
## mL.50year <- mean(sapply(regioni.50year,length)) #  11.2
## mH2.50year <- mean(test.50year["H2",]) #   1.53
## gH2gr2.50year <- sum(test.50year["H2",]>2)/696 #   0.34
## gH2gr4.50year <- sum(test.50year["H2",]>4)/696 #   0.07


###################################################
### code chunk number 13: nsRFA_ex02.Rnw:159-167 (eval = FALSE)
###################################################
## roi01.100year <- new.env()
## for(i in 1:696) {
##  print(paste(i,"/ 696"))
##  assign(as.character(row.names(roi.cd)[i]), roi.st.year(roi.cd[i,],as.data.frame(roi.cd),
##      row.names(roi.cd),am696[,"am"],am696[,"number"],test="HW",station.year=500,Nsim=100), 
##      env=roi01.100year)
## }
## roi01.100year <- as.list(roi01.100year)


###################################################
### code chunk number 14: nsRFA_ex02.Rnw:169-175 (eval = FALSE)
###################################################
## regioni.100year <- sapply(roi01.100year, estrai.region)
## test.100year <- sapply(roi01.100year, estrai.test)
## mL.100year <- mean(sapply(regioni.100year,length)) #  21.8
## mH2.100year <- mean(test.100year["H2",]) #   2.19
## gH2gr2.100year <- sum(test.100year["H2",]>2)/696 #   0.52
## gH2gr4.100year <- sum(test.100year["H2",]>4)/696 #   0.15


###################################################
### code chunk number 15: nsRFA_ex02.Rnw:178-184 (eval = FALSE)
###################################################
## table16.2 <- data.frame(signif(rbind(c(mL.50year,mH2.50year,
##                                        gH2gr2.50year*100,gH2gr4.50year*100),
##               c(mL.100year,mH2.100year,gH2gr2.100year*100,gH2gr4.100year*100)),3), 
##               row.names=c("50-year","100-year"))
## names(table16.2) <- c("Avg. n sites","m(H2)","% H2>2","% H2>4")
## print(table16.2)


###################################################
### code chunk number 16: nsRFA_ex02.Rnw:186-190
###################################################
table16.2 <- data.frame(signif(rbind(c(11.2,1.53,0.34*100,0.07*100),
              c(21.8,2.19,0.52*100,0.15*100)),3), row.names=c("50-year","100-year"))
names(table16.2) <- c("Avg. n sites","m(H2)","% H2>2","% H2>4")
print(table16.2)


###################################################
### code chunk number 17: nsRFA_ex02.Rnw:196-201
###################################################
prova54088 <- roi.st.year(roi.cd["54088",],roi.cd,row.names(roi.cd),am696[,"am"],
                         am696[,"number"],test="HW",station.year=250,Nsim=500)

prova28018 <- roi.st.year(roi.cd["28018",],roi.cd,row.names(roi.cd),am696[,"am"],
                          am696[,"number"],test="HW",station.year=250,Nsim=500)


###################################################
### code chunk number 18: nsRFA_ex02.Rnw:207-222
###################################################
Lmomenti696 <- as.data.frame(Lmomenti696)
par(mfrow=c(1,2))
 plot(Lmomenti696[c("lca","lcv")], xlab="L-CA", ylab="L-CV",
      pch=".", cex=2, main="54088"); grid()
 points(Lmomenti696[c("54088"), c("lca","lcv")],
        pch=19, col="red", cex=1)
 points(Lmomenti696[prova54088$region[-1], c("lca","lcv")],
        pch=19, cex=1)
 plot(Lmomenti696[,c("lca","lkur")], xlab="L-CA", ylab="L-kur",
      pch=".", cex=2, main="28018"); grid()
 points(Lmomenti696[c("28018"), c("lca","lcv")],
        pch=19, col="red", cex=1)
 points(Lmomenti696[prova28018$region[-1], c("lca","lcv")],
        pch=19, cex=1)
par(mfrow=c(1,1))


###################################################
### code chunk number 19: nsRFA_ex02.Rnw:236-296
###################################################
figure16.9a <- function (x,r,cd) {
 # x = station of interest (e.g. "28018")
 # r = output of roi.st.year()

 if(!r$region[1]==x) r$region <- c(x,r$region)
 row.names(cd) <- cd[,"number"]
 n <- length(cd[,"number"])
 cd.r <- cd[r$region,]
 par(mfrow=c(2,3))
  hist(log(cd[,"dtm_area"]),col="lightgray",border="lightgray",
       main="",xlab="AREA",axes=FALSE)
  axis(1,at=c(log(1),log(10),log(100),log(1000),log(10000)),
       label=c("1","10","100","1000","10000"))
  axis(2,at=seq(0,1,by=.05)*n,label=seq(0,1,by=.05))
  box()
  points(cbind(log(cd.r[-1,"dtm_area"]),0),pch=19,cex=.7)
  points(cbind(log(cd.r[1,"dtm_area"]),0),pch=4,cex=2,lwd=2)

  hist(cd[,"saar"],col="lightgray",border="lightgray",
       main="",xlab="SAAR",axes=FALSE)
  axis(1)
  axis(2,at=seq(0,1,by=.05)*n,label=seq(0,1,by=.05))
  box()
  points(cbind(cd.r[-1,"saar"],0),pch=19,cex=.7)
  points(cbind(cd.r[1,"saar"],0),pch=4,cex=2,lwd=2)

  hist(cd[,"bfihost"],col="lightgray",border="lightgray",
       main="",xlab="BFIHOST",axes=FALSE)
  axis(1)
  axis(2,at=seq(0,1,by=.05)*n,label=seq(0,1,by=.05))
  box()
  points(cbind(cd.r[-1,"bfihost"],0),pch=19,cex=.7)
  points(cbind(cd.r[1,"bfihost"],0),pch=4,cex=2,lwd=2)

  hist(cd[,"farl"],col="lightgray",border="lightgray",
       main="",xlab="FARL",axes=FALSE)
  axis(1)
  axis(2,at=seq(0,1,by=.05)*n,label=seq(0,1,by=.05))
  box()
  points(cbind(cd.r[-1,"farl"],0),pch=19,cex=.7)
  points(cbind(cd.r[1,"farl"],0),pch=4,cex=2,lwd=2)

  hist(cd[,"propwet"],col="lightgray",border="lightgray",
       main="",xlab="PROPWET",axes=FALSE)
  axis(1)
  axis(2,at=seq(0,1,by=.05)*n,label=seq(0,1,by=.05))
  box()
  points(cbind(cd.r[-1,"propwet"],0),pch=19,cex=.7)
  points(cbind(cd.r[1,"propwet"],0),pch=4,cex=2,lwd=2)

  hist(cd[,"urbext1990"],col="lightgray",border="lightgray",
       main="",xlab="URBEXT",axes=FALSE)
  axis(1)
  axis(2,at=seq(0,1,by=.05)*n,label=seq(0,1,by=.05))
  box()
  points(cbind(cd.r[-1,"urbext1990"],0),pch=19,cex=.7)
  points(cbind(cd.r[1,"urbext1990"],0),pch=4,cex=2,lwd=2)
 par(mfrow=c(1,1))
 title(main=x,cex.main=1,font.main=1)
}


###################################################
### code chunk number 20: nsRFA_ex02.Rnw:298-300
###################################################
prova40009 <- roi.st.year(roi.cd["40009",],roi.cd,row.names(roi.cd),am696[,"am"],
                          am696[,"number"],test="HW",station.year=500,Nsim=500)


###################################################
### code chunk number 21: nsRFA_ex02.Rnw:306-307
###################################################
figure16.9a("40009",prova40009,cd696)


###################################################
### code chunk number 22: nsRFA_ex02.Rnw:319-394
###################################################
figure16.9b <- function (x,r,am,cd) {
 # x = station of interest (e.g. "28018")
 # r = output of roi.st.year()

 row.names(cd) <- cd[,"number"]
 n <- length(cd[,"number"])
 cd.r <- cd[r$region,]
 cd.x <- cd[x,]
 fac <- factor(am[,"number"],levels=cd.r[,"number"])
 am.r <- am[!is.na(fac),]
 fac <- factor(am[,"number"],levels=x)
 am.x <- am[!is.na(fac),]
 am.xr <- rbind(am.x,am.r)
 QMED.r <- tapply(am.r[,4],am.r[,1],median)
 QMED.x <- median(am.x[,4])
 am.r.adim <- am.r; am.r.adim[,4] <- am.r[,4]/unsplit(QMED.r,am.r[,1])
 am.x.adim <- am.x; am.x.adim[,4] <- am.x[,4]/QMED.x
 lcv <- tapply(am[,4],am[,1],LCV)
 lca <- tapply(am[,4],am[,1],LCA)
 lkur <- tapply(am[,4],am[,1],Lkur)
 lcv.r <- tapply(am.r[,4],am.r[,1],LCV)
 lca.r <- tapply(am.r[,4],am.r[,1],LCA)
 lkur.r <- tapply(am.r[,4],am.r[,1],Lkur)
 lcv.x <- LCV(am.x[,4])
 lca.x <- LCA(am.x[,4])
 lkur.x <- Lkur(am.x[,4])
 days <- as.numeric(format(as.Date(am[,2]),"%j"))
 days.r <- as.numeric(format(as.Date(am.r[,2]),"%j"))
 days.x <- as.numeric(format(as.Date(am.x[,2]),"%j"))

 par(mfrow=c(2,3))
  lognormplot(am.r.adim[,4],line=FALSE,xlab="Q/QMED",type="n")
  for(i in r$region) {
   xxx <- am.r.adim[am.r.adim[,1]==i,4]
   normpoints(xxx,type="l",col="gray")
  }
  normpoints(am.r.adim[,4],type="l",lwd=2)
  normpoints(am.x.adim[,4],type="l",col=2,lwd=2)

  plot(lca,lcv,pch=".",cex=2)
  points(lca.r,lcv.r,pch=19)
  points(lca.x,lcv.x,pch=4,cex=2,lwd=2)

  plot(lca,lkur,pch=".",cex=2)
  points(lca.r,lkur.r,pch=19)
  points(lca.x,lkur.x,pch=4,cex=2,lwd=2)

  plot(cd[c("ihdtm_ngr_x","ihdtm_ngr_y")],pch=".",cex=2,xlab="",ylab="",axes=FALSE)
  points(cd.r[c("ihdtm_ngr_x","ihdtm_ngr_y")],pch=19)
  points(cd.x[c("ihdtm_ngr_x","ihdtm_ngr_y")],pch=4,cex=2,lwd=2)

  consistencyplot (am.r[,3],am.r[,1])

  dummy <- seq(0,2*pi,length=100)
  plot(cos(dummy),sin(dummy),type="l",xlab="",ylab="",axes=FALSE)
  abline(h=0,lty=3); abline(v=0,lty=3)
  radd <- days*pi/180
  XFLOOD <- tapply(cos(radd),am[,1],mean)
  YFLOOD <- tapply(sin(radd),am[,1],mean)
  points(XFLOOD,YFLOOD,pch=".",cex=2)
  radd <- days.r*pi/180
  XFLOOD <- tapply(cos(radd),am.r[,1],mean)
  YFLOOD <- tapply(sin(radd),am.r[,1],mean)
  points(XFLOOD,YFLOOD,pch=19,cex=1)
  radd <- days.x*pi/180
  XFLOOD <- tapply(cos(radd),am.x[,1],mean)
  YFLOOD <- tapply(sin(radd),am.x[,1],mean)
  points(XFLOOD,YFLOOD,pch=4,cex=2,lwd=2)
  axis(1,at=0,label="Oct 1")
  axis(2,at=0,label="Jul 1")
  axis(3,at=0,label="Apr 1")
  axis(4,at=0,label="Jan 1")
 par(mfrow=c(1,1))
 title(main=x,cex.main=1,font.main=1)
}


###################################################
### code chunk number 23: nsRFA_ex02.Rnw:400-401
###################################################
figure16.9b("40009",prova40009,am696,cd696)


###################################################
### code chunk number 24: nsRFA_ex02.Rnw:418-420
###################################################
prova45001 <- roi.st.year(roi.cd["45001",],roi.cd,row.names(roi.cd),am696[,"am"],
                          am696[,"number"],test="HW",station.year=250,Nsim=500)


###################################################
### code chunk number 25: nsRFA_ex02.Rnw:426-427
###################################################
figure16.9a("45001",prova45001,cd696)


###################################################
### code chunk number 26: nsRFA_ex02.Rnw:442-443
###################################################
figure16.9b("45001",prova45001,am696,cd696)


