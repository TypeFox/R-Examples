### R code from vignette source 'austides.Rnw'

###################################################
### code chunk number 1: austides.Rnw:108-112
###################################################
years <- 2012:2014
sites <- c("CapeFerguson","PortKembla","Portland","Thevenard",
           "Esperance","Hillarys","Broome","Darwin")
siteid <- paste0("IDO710",c("01","03","08","10","11","12","13","14"))


###################################################
### code chunk number 2: austides.Rnw:115-116
###################################################
options(width=60)


###################################################
### code chunk number 3: austides.Rnw:121-137 (eval = FALSE)
###################################################
## fn <- "http://www.bom.gov.au/ntc/"
## fn <- paste0(fn, siteid,"/",siteid,"_")
## abslmp <- list(); tp <- tempfile()
## for(i in 1:length(sites)) {
##   nms <- paste0(fn[i], years,".csv")
##   abslmp[[i]] <- data.frame()
##   for(j in 1:length(years)) {
##     download.file(nms[j], tp)
##     abslmp[[i]] <- rbind(abslmp[[i]], read.csv(tp,as.is=TRUE)[,1:2])
##   }
##   colnames(abslmp[[i]]) <- c("DateTime","SeaLevel")
##   abslmp[[i]]$DateTime <- as.POSIXct(abslmp[[i]]$DateTime, tz="UTC", 
##     format = "%d-%b-%Y %H:%M")
##   abslmp[[i]]$SeaLevel[abslmp[[i]]$SeaLevel == -9999] <- NA
## }
## names(abslmp) <- sites


###################################################
### code chunk number 4: austides.Rnw:144-148
###################################################
library(TideHarmonics)
abslmp <- list()
for(i in 1:length(sites)) abslmp[[i]] <- get(sites[i])
names(abslmp) <- sites


###################################################
### code chunk number 5: austides.Rnw:157-161
###################################################
library(TideHarmonics)
esp <- abslmp$Esperance
m1 <- ftide(esp$SeaLevel, esp$DateTime, hc4)
m1


###################################################
### code chunk number 6: austides.Rnw:194-198 (eval = FALSE)
###################################################
## t1 <- as.POSIXct("2014-12-01 00:00", tz = "UTC")
## t2 <- as.POSIXct("2014-12-31 00:00", tz = "UTC")
## plot(m1, t1, t2)
## plot(m1, t1, t2, split = TRUE, ylim = c(-0.18,0.18))


###################################################
### code chunk number 7: esp
###################################################
t1 <- as.POSIXct("2014-12-01 00:00", tz = "UTC")
t2 <- as.POSIXct("2014-12-31 00:00", tz = "UTC")
plot(m1, t1, t2, cex.lab = 1.5, cex.axis = 1.5)


###################################################
### code chunk number 8: espK1
###################################################
t1 <- as.POSIXct("2014-12-01 00:00", tz = "UTC")
t2 <- as.POSIXct("2014-12-31 00:00", tz = "UTC")
plot(m1, t1, t2, split = TRUE, ylim = c(-0.18,0.18), which = "K1", cex.lab = 1.5, cex.axis = 1.5,cex.main = 1.25)


###################################################
### code chunk number 9: espS2
###################################################
t1 <- as.POSIXct("2014-12-01 00:00", tz = "UTC")
t2 <- as.POSIXct("2014-12-31 00:00", tz = "UTC")
plot(m1, t1, t2, split = TRUE, ylim = c(-0.18,0.18), which = "S2", cex.lab = 1.5, cex.axis = 1.5,cex.main = 1.25)


###################################################
### code chunk number 10: espO1
###################################################
t1 <- as.POSIXct("2014-12-01 00:00", tz = "UTC")
t2 <- as.POSIXct("2014-12-31 00:00", tz = "UTC")
plot(m1, t1, t2, split = TRUE, ylim = c(-0.18,0.18), which = "O1", cex.lab = 1.5, cex.axis = 1.5,cex.main = 1.25)


###################################################
### code chunk number 11: espM2
###################################################
t1 <- as.POSIXct("2014-12-01 00:00", tz = "UTC")
t2 <- as.POSIXct("2014-12-31 00:00", tz = "UTC")
plot(m1, t1, t2, split = TRUE, ylim = c(-0.18,0.18), which = "M2", cex.lab = 1.5, cex.axis = 1.5,cex.main = 1.25)


###################################################
### code chunk number 12: austides.Rnw:233-234
###################################################
t1 <- as.POSIXct("2014-12-01 00:00", tz = "UTC")
t2 <- as.POSIXct("2014-12-31 00:00", tz = "UTC")
plot(m1, t1, t2, cex.lab = 1.5, cex.axis = 1.5)


###################################################
### code chunk number 13: austides.Rnw:244-245
###################################################
t1 <- as.POSIXct("2014-12-01 00:00", tz = "UTC")
t2 <- as.POSIXct("2014-12-31 00:00", tz = "UTC")
plot(m1, t1, t2, split = TRUE, ylim = c(-0.18,0.18), which = "K1", cex.lab = 1.5, cex.axis = 1.5,cex.main = 1.25)


###################################################
### code chunk number 14: austides.Rnw:247-248
###################################################
t1 <- as.POSIXct("2014-12-01 00:00", tz = "UTC")
t2 <- as.POSIXct("2014-12-31 00:00", tz = "UTC")
plot(m1, t1, t2, split = TRUE, ylim = c(-0.18,0.18), which = "S2", cex.lab = 1.5, cex.axis = 1.5,cex.main = 1.25)


###################################################
### code chunk number 15: austides.Rnw:250-251
###################################################
t1 <- as.POSIXct("2014-12-01 00:00", tz = "UTC")
t2 <- as.POSIXct("2014-12-31 00:00", tz = "UTC")
plot(m1, t1, t2, split = TRUE, ylim = c(-0.18,0.18), which = "O1", cex.lab = 1.5, cex.axis = 1.5,cex.main = 1.25)


###################################################
### code chunk number 16: austides.Rnw:253-254
###################################################
t1 <- as.POSIXct("2014-12-01 00:00", tz = "UTC")
t2 <- as.POSIXct("2014-12-31 00:00", tz = "UTC")
plot(m1, t1, t2, split = TRUE, ylim = c(-0.18,0.18), which = "M2", cex.lab = 1.5, cex.axis = 1.5,cex.main = 1.25)


###################################################
### code chunk number 17: austides.Rnw:272-281
###################################################
m2 <- ftide(esp$SeaLevel, esp$DateTime, hc60)
tt <- c(rep(2014,12),2015)
tt <- paste(tt,sprintf("%02d",c(1:12,1)),"01",sep="-")
tt <- as.POSIXct(tt, tz = "UTC")
for(i in 1:12) {
  plot(m2, tt[i], tt[i+1], main=paste(month.abb[i],2014))
  ind <- esp$DateTime >= tt[i] & esp$DateTime <= tt[i+1]
  lines(esp[ind,], lty=2, col="red")
}


###################################################
### code chunk number 18: esp1
###################################################
i <- 1
plot(m2, tt[i], tt[i+1], main=paste(month.abb[i],2014),cex.lab = 1.5, cex.axis = 1.5,cex.main = 1.25)
ind <- esp$DateTime >= tt[i] & esp$DateTime <= tt[i+1]
lines(esp[ind,], lty=2, col="red")


###################################################
### code chunk number 19: esp2
###################################################
i <- 2
plot(m2, tt[i], tt[i+1], main=paste(month.abb[i],2014),cex.lab = 1.5, cex.axis = 1.5,cex.main = 1.25)
ind <- esp$DateTime >= tt[i] & esp$DateTime <= tt[i+1]
lines(esp[ind,], lty=2, col="red")


###################################################
### code chunk number 20: esp3
###################################################
i <- 3
plot(m2, tt[i], tt[i+1], main=paste(month.abb[i],2014),cex.lab = 1.5,cex.axis = 1.5,cex.main = 1.25)
ind <- esp$DateTime >= tt[i] & esp$DateTime <= tt[i+1]
lines(esp[ind,], lty=2, col="red")


###################################################
### code chunk number 21: esp4
###################################################
i <- 4
plot(m2, tt[i], tt[i+1], main=paste(month.abb[i],2014),cex.lab = 1.5,cex.axis = 1.5,cex.main = 1.25)
ind <- esp$DateTime >= tt[i] & esp$DateTime <= tt[i+1]
lines(esp[ind,], lty=2, col="red")


###################################################
### code chunk number 22: austides.Rnw:314-315
###################################################
i <- 1
plot(m2, tt[i], tt[i+1], main=paste(month.abb[i],2014),cex.lab = 1.5, cex.axis = 1.5,cex.main = 1.25)
ind <- esp$DateTime >= tt[i] & esp$DateTime <= tt[i+1]
lines(esp[ind,], lty=2, col="red")


###################################################
### code chunk number 23: austides.Rnw:317-318
###################################################
i <- 2
plot(m2, tt[i], tt[i+1], main=paste(month.abb[i],2014),cex.lab = 1.5, cex.axis = 1.5,cex.main = 1.25)
ind <- esp$DateTime >= tt[i] & esp$DateTime <= tt[i+1]
lines(esp[ind,], lty=2, col="red")


###################################################
### code chunk number 24: austides.Rnw:320-321
###################################################
i <- 3
plot(m2, tt[i], tt[i+1], main=paste(month.abb[i],2014),cex.lab = 1.5,cex.axis = 1.5,cex.main = 1.25)
ind <- esp$DateTime >= tt[i] & esp$DateTime <= tt[i+1]
lines(esp[ind,], lty=2, col="red")


###################################################
### code chunk number 25: austides.Rnw:323-324
###################################################
i <- 4
plot(m2, tt[i], tt[i+1], main=paste(month.abb[i],2014),cex.lab = 1.5,cex.axis = 1.5,cex.main = 1.25)
ind <- esp$DateTime >= tt[i] & esp$DateTime <= tt[i+1]
lines(esp[ind,], lty=2, col="red")


###################################################
### code chunk number 26: austides.Rnw:344-352 (eval = FALSE)
###################################################
## mlst <- list()
## for(i in 1:length(sites)) {
##   df <- abslmp[[i]]
##   mlst[[i]] <- ftide(df$SeaLevel, df$DateTime)
## }
## names(mlst) <- sites
## sapply(mlst, function(x) x[["fval"]])
## lapply(mlst, function(x) round(head(coef(x, hc = TRUE),10),3))


###################################################
### code chunk number 27: austides.Rnw:355-362
###################################################
mlst <- list()
for(i in 1:length(sites)) {
  df <- abslmp[[i]]
  mlst[[i]] <- ftide(df$SeaLevel, df$DateTime)
}
names(mlst) <- sites
sapply(mlst, function(x) x[["fval"]])


###################################################
### code chunk number 28: austides.Rnw:371-380
###################################################
t1 <- as.POSIXct("2014-12-01 00:00", tz = "UTC")
t2 <- as.POSIXct("2014-12-31 00:00", tz = "UTC")
msl <- sapply(mlst, function(x) x[["msl"]])
for(i in 1:length(sites)) {
  plot(mlst[[i]], t1, t2, msl = FALSE, main=sites[i])
  df <- abslmp[[i]]
  ind <- df$DateTime >= t1 & df$DateTime <= t2
  lines(df[ind,1], df[ind,2]-msl[i], lty=2, col="red")
}


###################################################
### code chunk number 29: site1
###################################################
i <- 1
plot(mlst[[i]], t1, t2, msl = FALSE, main=sites[i],cex.lab = 1.5, cex.axis = 1.5,cex.main = 1.25)
df <- abslmp[[i]]
ind <- df$DateTime >= t1 & df$DateTime <= t2
lines(df[ind,1], df[ind,2]-msl[i], lty=2, col="red")


###################################################
### code chunk number 30: site2
###################################################
i <- 2
plot(mlst[[i]], t1, t2, msl = FALSE, main=sites[i],cex.lab = 1.5, cex.axis = 1.5,cex.main = 1.25)
df <- abslmp[[i]]
ind <- df$DateTime >= t1 & df$DateTime <= t2
lines(df[ind,1], df[ind,2]-msl[i], lty=2, col="red")


###################################################
### code chunk number 31: site3
###################################################
i <- 3
plot(mlst[[i]], t1, t2, msl = FALSE, main=sites[i],cex.lab = 1.5, cex.axis = 1.5,cex.main = 1.25)
df <- abslmp[[i]]
ind <- df$DateTime >= t1 & df$DateTime <= t2
lines(df[ind,1], df[ind,2]-msl[i], lty=2, col="red")


###################################################
### code chunk number 32: site4
###################################################
i <- 4
plot(mlst[[i]], t1, t2, msl = FALSE, main=sites[i],cex.lab = 1.5, cex.axis = 1.5,cex.main = 1.25)
df <- abslmp[[i]]
ind <- df$DateTime >= t1 & df$DateTime <= t2
lines(df[ind,1], df[ind,2]-msl[i], lty=2, col="red")


###################################################
### code chunk number 33: site5
###################################################
i <- 5
plot(mlst[[i]], t1, t2, msl = FALSE, main=sites[i],cex.lab = 1.5, cex.axis = 1.5,cex.main = 1.25)
df <- abslmp[[i]]
ind <- df$DateTime >= t1 & df$DateTime <= t2
lines(df[ind,1], df[ind,2]-msl[i], lty=2, col="red")


###################################################
### code chunk number 34: site6
###################################################
i <- 6
plot(mlst[[i]], t1, t2, msl = FALSE, main=sites[i],cex.lab = 1.5, cex.axis = 1.5,cex.main = 1.25)
df <- abslmp[[i]]
ind <- df$DateTime >= t1 & df$DateTime <= t2
lines(df[ind,1], df[ind,2]-msl[i], lty=2, col="red")


###################################################
### code chunk number 35: site7
###################################################
i <- 7
plot(mlst[[i]], t1, t2, msl = FALSE, main=sites[i],cex.lab = 1.5, cex.axis = 1.5,cex.main = 1.25)
df <- abslmp[[i]]
ind <- df$DateTime >= t1 & df$DateTime <= t2
lines(df[ind,1], df[ind,2]-msl[i], lty=2, col="red")


###################################################
### code chunk number 36: site8
###################################################
i <- 8
plot(mlst[[i]], t1, t2, msl = FALSE, main=sites[i],cex.lab = 1.5, cex.axis = 1.5,cex.main = 1.25)
df <- abslmp[[i]]
ind <- df$DateTime >= t1 & df$DateTime <= t2
lines(df[ind,1], df[ind,2]-msl[i], lty=2, col="red")


###################################################
### code chunk number 37: austides.Rnw:449-450
###################################################
i <- 1
plot(mlst[[i]], t1, t2, msl = FALSE, main=sites[i],cex.lab = 1.5, cex.axis = 1.5,cex.main = 1.25)
df <- abslmp[[i]]
ind <- df$DateTime >= t1 & df$DateTime <= t2
lines(df[ind,1], df[ind,2]-msl[i], lty=2, col="red")


###################################################
### code chunk number 38: austides.Rnw:452-453
###################################################
i <- 2
plot(mlst[[i]], t1, t2, msl = FALSE, main=sites[i],cex.lab = 1.5, cex.axis = 1.5,cex.main = 1.25)
df <- abslmp[[i]]
ind <- df$DateTime >= t1 & df$DateTime <= t2
lines(df[ind,1], df[ind,2]-msl[i], lty=2, col="red")


###################################################
### code chunk number 39: austides.Rnw:455-456
###################################################
i <- 3
plot(mlst[[i]], t1, t2, msl = FALSE, main=sites[i],cex.lab = 1.5, cex.axis = 1.5,cex.main = 1.25)
df <- abslmp[[i]]
ind <- df$DateTime >= t1 & df$DateTime <= t2
lines(df[ind,1], df[ind,2]-msl[i], lty=2, col="red")


###################################################
### code chunk number 40: austides.Rnw:458-459
###################################################
i <- 4
plot(mlst[[i]], t1, t2, msl = FALSE, main=sites[i],cex.lab = 1.5, cex.axis = 1.5,cex.main = 1.25)
df <- abslmp[[i]]
ind <- df$DateTime >= t1 & df$DateTime <= t2
lines(df[ind,1], df[ind,2]-msl[i], lty=2, col="red")


###################################################
### code chunk number 41: austides.Rnw:469-470
###################################################
i <- 5
plot(mlst[[i]], t1, t2, msl = FALSE, main=sites[i],cex.lab = 1.5, cex.axis = 1.5,cex.main = 1.25)
df <- abslmp[[i]]
ind <- df$DateTime >= t1 & df$DateTime <= t2
lines(df[ind,1], df[ind,2]-msl[i], lty=2, col="red")


###################################################
### code chunk number 42: austides.Rnw:472-473
###################################################
i <- 6
plot(mlst[[i]], t1, t2, msl = FALSE, main=sites[i],cex.lab = 1.5, cex.axis = 1.5,cex.main = 1.25)
df <- abslmp[[i]]
ind <- df$DateTime >= t1 & df$DateTime <= t2
lines(df[ind,1], df[ind,2]-msl[i], lty=2, col="red")


###################################################
### code chunk number 43: austides.Rnw:475-476
###################################################
i <- 7
plot(mlst[[i]], t1, t2, msl = FALSE, main=sites[i],cex.lab = 1.5, cex.axis = 1.5,cex.main = 1.25)
df <- abslmp[[i]]
ind <- df$DateTime >= t1 & df$DateTime <= t2
lines(df[ind,1], df[ind,2]-msl[i], lty=2, col="red")


###################################################
### code chunk number 44: austides.Rnw:478-479
###################################################
i <- 8
plot(mlst[[i]], t1, t2, msl = FALSE, main=sites[i],cex.lab = 1.5, cex.axis = 1.5,cex.main = 1.25)
df <- abslmp[[i]]
ind <- df$DateTime >= t1 & df$DateTime <= t2
lines(df[ind,1], df[ind,2]-msl[i], lty=2, col="red")


