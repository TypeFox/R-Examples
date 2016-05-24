### R code from vignette source 'simulations.rnw'

###################################################
### code chunk number 1: Setup
###################################################

rm(list=ls())

options(width=65)

if(!require(lattice, quietly=TRUE)) install.packages("lattice")
if(!require(Hmisc, quietly=TRUE)) install.packages("Hmisc")
if(!require(rconifers, quietly=TRUE)) install.packages("rconifers")

## this is for
b201.age <- 2:15*10
b201.tht <- c(47,83,108,125,140,152,163,172,180,187,192,196,199,201)
b201.tpa <- c(756,483,335,248,195,160,136,118,106,95,87,80,75,71)
b201.qmd <- c(4.9,7.6,10.2,12.8,15.2,17.5,19.6,21.4,23.1,24.6,26.1,27.5,28.8,30)
b201.ba <- c(99,152,191,220,244,264,280,294,306,317,326,335,343,351)
b201.vol <- c(1650,4330,6900,9320,11450,13300,14990,16400,17550,18510,19320,20000,20640,21270)
b201.mai <- c(77,136,164,177,181,181,178,173) ## only goes from 20 to 90 age[1:8]



###################################################
### code chunk number 2: simulations.rnw:213-214 (eval = FALSE)
###################################################
## install.packages("rconifers")


###################################################
### code chunk number 3: simulations.rnw:220-224
###################################################

library(rconifers)
set.species.map(set.variant(1))



###################################################
### code chunk number 4: Maybe should be ... (eval = FALSE)
###################################################
## set.variant(1)
## data( smc )
## set.species.map( smc )


###################################################
### code chunk number 5: simulations.rnw:238-241
###################################################
dim(smc)
names(smc)
head(smc)


###################################################
### code chunk number 6: simulations.rnw:281-285
###################################################

data(plots.smc)
data(plants.smc)



###################################################
### code chunk number 7: simulations.rnw:291-298
###################################################

sample.3 <- list(plots = plots.smc, 
                 plants = plants.smc, 
                 age = 3, 
                 x0 = 0.0)
class(sample.3)  <- "sample.data"



###################################################
### code chunk number 8: simulations.rnw:304-307
###################################################

summary(sample.3)



###################################################
### code chunk number 9: simulations.rnw:325-333
###################################################
summary(sample.3)
sample.25 <- project(sample.3, 
                     22,
                     control = list(rand.err = 0,
                                    rand.seed = 0,
                                    endemic.mort = 0,
                                    sdi.mort = 0))
summary(sample.25)


###################################################
### code chunk number 10: simulations.rnw:371-391
###################################################

res.v <- vector(length = 98, mode = "list")
s0 <- sample.3
res.v[[1]] <- data.frame(age = s0$age, sp.sums(s0)["DF",])

## grow from age 3 to age 100 (97 more years)
for(m in 1:97) {
  ## project s0 in one year intervals
  s1 <- project(s0, 1, control = list(rand.err = 0,
                                      rand.seed = 0,
                                      endemic.mort = 1,
                                      sdi.mort = 1))          
  res.v[[m+1]] <- data.frame(age = s1$age, 
                             sp.sums(s1)["DF",])
  s0 <- s1      
}

res.v <- do.call(rbind, res.v)
res <- res.v



###################################################
### code chunk number 11: fig-rconifers
###################################################

opar <- par(mfcol=c(1,2), las=1, cex=0.8)

################################################################################
## plot the trees per acre, 7.0 inches dbh and larger
## table three (3). 
## reflects pure even-aged second growth stand
##plot(rcon.111$tpa.7.plus ~ rcon.111$age, 
#plot(rcon.111$expf ~ rcon.111$age, 
plot(res$expf ~ res$age, 
     type="n", 
     ylim=c(0,350),
     xlim=c(0,100),
##      main="Trees per Acre > 7 Inches",
##      ylab="Trees/Acre > 7.0 Inches", 
      main="Trees per Acre (All Stems)",
      ylab="Trees/Acre (All Stems)", 
     xlab="Total Age, in years", col="black")
grid(lty=1)
par(new=TRUE)
##plot(rcon.111$tpa.7.plus ~ rcon.111$age, 
##plot(rcon.111$expf ~ rcon.111$age, 
plot(res$expf ~ res$age, 
     type="l", 
     ylim=c(0,350),
     xlim=c(0,100),
     main=NA,
     ylab=NA, 
     xlab=NA, col="black")
par(new=TRUE)
##plot(rcon.531$tpa.7.plus ~ rcon.531$age, 
##plot(rcon.531$expf ~ rcon.531$age, 
plot(res$expf ~ res$age, 
     type="l", 
     ylim=c(0,350),
     xlim=c(0,100),
     main=NA,
     ylab=NA, 
     xlab=NA, col="black")
lines(b201.age, b201.tpa, lty=2, col="gray")


################################################################################
## plot the basal area 
##plot(rcon.111$ba ~ rcon.111$age, 
plot(res$ba ~ res$age, 
     type="n", 
     ylim=c(0,300),
     xlim=c(0,100),
     main="Basal Area",
     ylab="Ft^2/Acre", 
     xlab="Total Age, in years", col="black")
grid(lty=1)
par(new=TRUE)
##plot(rcon.111$ba ~ rcon.111$age, 
plot(res$ba ~ res$age, 
     type="l", 
     ylim=c(0,300),
     xlim=c(0,100),
     main="Basal Area",
     ylab="Ft^2/Acre", 
     xlab="Total Age, in years", col="black")
par(new=TRUE)
##plot(rcon.531$ba ~ rcon.531$age, 
plot(res$ba ~ res$age, 
     type="l", 
     ylim=c(0,300),
     xlim=c(0,100),
     main=NA,
     ylab=NA, 
     xlab=NA, col="black")
lines(b201.age, b201.ba, lty=2, col="gray")



###################################################
### code chunk number 12: fig-rconifers
###################################################

opar <- par(mfcol=c(1,2), las=1, cex=0.8)

################################################################################
## plot the trees per acre, 7.0 inches dbh and larger
## table three (3). 
## reflects pure even-aged second growth stand
##plot(rcon.111$tpa.7.plus ~ rcon.111$age, 
#plot(rcon.111$expf ~ rcon.111$age, 
plot(res$expf ~ res$age, 
     type="n", 
     ylim=c(0,350),
     xlim=c(0,100),
##      main="Trees per Acre > 7 Inches",
##      ylab="Trees/Acre > 7.0 Inches", 
      main="Trees per Acre (All Stems)",
      ylab="Trees/Acre (All Stems)", 
     xlab="Total Age, in years", col="black")
grid(lty=1)
par(new=TRUE)
##plot(rcon.111$tpa.7.plus ~ rcon.111$age, 
##plot(rcon.111$expf ~ rcon.111$age, 
plot(res$expf ~ res$age, 
     type="l", 
     ylim=c(0,350),
     xlim=c(0,100),
     main=NA,
     ylab=NA, 
     xlab=NA, col="black")
par(new=TRUE)
##plot(rcon.531$tpa.7.plus ~ rcon.531$age, 
##plot(rcon.531$expf ~ rcon.531$age, 
plot(res$expf ~ res$age, 
     type="l", 
     ylim=c(0,350),
     xlim=c(0,100),
     main=NA,
     ylab=NA, 
     xlab=NA, col="black")
lines(b201.age, b201.tpa, lty=2, col="gray")


################################################################################
## plot the basal area 
##plot(rcon.111$ba ~ rcon.111$age, 
plot(res$ba ~ res$age, 
     type="n", 
     ylim=c(0,300),
     xlim=c(0,100),
     main="Basal Area",
     ylab="Ft^2/Acre", 
     xlab="Total Age, in years", col="black")
grid(lty=1)
par(new=TRUE)
##plot(rcon.111$ba ~ rcon.111$age, 
plot(res$ba ~ res$age, 
     type="l", 
     ylim=c(0,300),
     xlim=c(0,100),
     main="Basal Area",
     ylab="Ft^2/Acre", 
     xlab="Total Age, in years", col="black")
par(new=TRUE)
##plot(rcon.531$ba ~ rcon.531$age, 
plot(res$ba ~ res$age, 
     type="l", 
     ylim=c(0,300),
     xlim=c(0,100),
     main=NA,
     ylab=NA, 
     xlab=NA, col="black")
lines(b201.age, b201.ba, lty=2, col="gray")



###################################################
### code chunk number 13: shlib
###################################################

system("R CMD SHLIB chambers-1980.c")



###################################################
### code chunk number 14: simulations.rnw:703-706
###################################################

dyn.load("chambers-1980.so")



###################################################
### code chunk number 15: simulations.rnw:727-736
###################################################

site.index <- 120.0
total.age <- 60

.C("chambers_1980_stand_mean_height",
              as.double(site.index),
              as.double(total.age),
              smh = as.double(0))$smh



###################################################
### code chunk number 16: simulations.rnw:748-757
###################################################

chambers.1980.smh <- function(site.index,total.age) {
  ret.val <- .C("chambers_1980_stand_mean_height",
                as.double(site.index),
                as.double(total.age),
                smh = as.double(0))$smh
  ret.val
}



###################################################
### code chunk number 17: simulations.rnw:766-769
###################################################

chambers.1980.smh(120.0, 60.0)



###################################################
### code chunk number 18: simulations.rnw:782-784
###################################################
source("chambers-1980.r")



###################################################
### code chunk number 19: simulations.rnw:802-821
###################################################

chambers.1980 <- function(ages = 1:100, 
                          site = 125.0, 
                          pnba = 1.0) {
  ret.val <- matrix(0, length(ages), 5)
  for(i in ages) {    
    res <- c(i, 
             chambers.1980.adbh(site, i, pnba), ## Table 5
             chambers.1980.smh(site, i),        ## Table 16
             chambers.1980.nba(site, i),        ## Table 1
             chambers.1980.ntpa(site, i, pnba)  ## Table 3
            )
    ret.val[i,] <- res
  }
  ret.val <- as.data.frame(ret.val)
  names(ret.val) <- c("age","qmd","tht","ba","expf")
  ret.val
}



###################################################
### code chunk number 20: fig-chambers-1980
###################################################

## you need to convert between kings and jims

ch.site <- 125.0
pnba <- 1.0
ch.111 <- chambers.1980(ages=1:100, site=ch.site, pnba=pnba)
ch.531 <- chambers.1980(ages=1:100, site=ch.site, pnba=pnba)


## this code is not presented
opar <- par(mfcol=c(1,2), las=1, cex=0.8)
##opar <- par(mfcol=c(2,2), las=1, cex=0.8)


################################################################################
## plot the trees per acre, 7.0 inches dbh and larger
## table three (3). 
## reflects pure even-aged second growth stand
plot(ch.111$expf ~ ch.111$age, 
     type="n", 
     ylim=c(0,350),
     xlim=c(0,100),
     main="Trees per Acre",
     ylab="Trees/Acre > 7.0 Inches", 
     xlab="Total Age, in years")
grid(lty=1)
par(new=TRUE)
plot(ch.111$expf ~ ch.111$age, 
     type="l", col="darkgray",
     ylim=c(0,350),
     xlim=c(0,100),
     main=NA,
     ylab=NA, 
     xlab=NA)
## par(new=TRUE)
## plot(ch.531$expf ~ ch.531$age, 
##      type="l", col="darkgray",
##      ylim=c(0,350),
##      xlim=c(0,100),
##      main="Trees per Acre",
##      ylab="Trees/Acre > 7.0 Inches", 
##      xlab="Total Age, in years")
################################################################################
## plot the trees per acre, 7.0 inches dbh and larger
## table three (3). 
## reflects pure even-aged second growth stand
par(new=TRUE)
##plot(rcon.111$tpa.7.plus ~ rcon.111$age, 
##plot(rcon.111$expf ~ rcon.111$age, 
plot(res$expf ~ res$age, 
     type="l", col="black",
     ylim=c(0,350),
     xlim=c(0,100),
     main=NA,
     ylab=NA, 
     xlab=NA)
## par(new=TRUE)
## ##plot(rcon.531$tpa.7.plus ~ rcon.531$age, 
## plot(rcon.531$expf ~ rcon.531$age, 
##      type="l", col="black",
##      ylim=c(0,350),
##      xlim=c(0,100),
##      main=NA,
##      ylab=NA, 
##      xlab=NA)
lines(b201.age, b201.tpa, lty=2, col="darkgray")

################################################################################
## plot the basal area for the two models
################################################################################
plot(ch.111$ba ~ ch.111$age, 
     type="n", 
     ylim=c(0,300),
     xlim=c(0,100),
     col="black",
     main="Basal Area",
     ylab="Ft^2/Acre", 
     xlab="Total Age, in years")
grid(lty=1)
par(new=TRUE)
plot(ch.111$ba ~ ch.111$age, 
     type="l", 
     ylim=c(0,300),
     xlim=c(0,100),
     col="darkgray",
     main=NA,
     ylab=NA, 
     xlab=NA)
par(new=TRUE)
## plot(ch.531$ba ~ ch.531$age, 
##      type="l", 
##      ylim=c(0,300),
##      xlim=c(0,100),
##      col="darkgray",
##      main=NA,
##      ylab=NA, 
##      xlab=NA)
## par(new=TRUE)
##plot(rcon.111$ba ~ rcon.111$age, 
plot(res$ba ~ res$age, 
     type="l", 
     ylim=c(0,300),
     xlim=c(0,100),
     col="black",
     main=NA,
     ylab=NA, 
     xlab=NA)
## par(new=TRUE)
## plot(rcon.531$ba ~ rcon.531$age, 
##      type="l", 
##      ylim=c(0,300),
##      xlim=c(0,100),
##      col="black",
##      main=NA,
##      ylab=NA, 
##      xlab=NA)
lines(b201.age, b201.ba, lty=2, col="darkgray")


################################################################################
## plot the quadratic mean diameter for the two models
################################################################################
## plot(ch.111$qmd ~ ch.111$age, 
##      type="n", 
##      ylim=c(0,25),
##      xlim=c(0,100),
##      col="black",
##      main="Mean Stem Diameter",
##      ylab="Inches", 
##      xlab="Total Age, in years")
## grid(lty=1)
## par(new=TRUE)
## plot(ch.111$qmd ~ ch.111$age, 
##      type="l", 
##      ylim=c(0,25),
##      xlim=c(0,100),
##      col="darkgray",
##      main=NA,
##      ylab=NA, 
##      xlab=NA)
## par(new=TRUE)
## plot(ch.531$qmd ~ ch.531$age, 
##      type="l", 
##      ylim=c(0,25),
##      xlim=c(0,100),
##      col="darkgray",
##      main=NA,
##      ylab=NA, 
##      xlab=NA)
## par(new=TRUE)
## plot(rcon.111$qmd ~ rcon.111$age, 
##      type="l", 
##      ylim=c(0,25),
##      xlim=c(0,100),
##      col="black",
##      main=NA,
##      ylab=NA, 
##      xlab=NA)
## par(new=TRUE)
## plot(rcon.531$qmd ~ rcon.531$age, 
##      type="l", 
##      ylim=c(0,25),
##      xlim=c(0,100),
##      col="black",
##      main=NA,
##      ylab=NA, 
##      xlab=NA)



###################################################
### code chunk number 21: fig-chambers-1980
###################################################

## you need to convert between kings and jims

ch.site <- 125.0
pnba <- 1.0
ch.111 <- chambers.1980(ages=1:100, site=ch.site, pnba=pnba)
ch.531 <- chambers.1980(ages=1:100, site=ch.site, pnba=pnba)


## this code is not presented
opar <- par(mfcol=c(1,2), las=1, cex=0.8)
##opar <- par(mfcol=c(2,2), las=1, cex=0.8)


################################################################################
## plot the trees per acre, 7.0 inches dbh and larger
## table three (3). 
## reflects pure even-aged second growth stand
plot(ch.111$expf ~ ch.111$age, 
     type="n", 
     ylim=c(0,350),
     xlim=c(0,100),
     main="Trees per Acre",
     ylab="Trees/Acre > 7.0 Inches", 
     xlab="Total Age, in years")
grid(lty=1)
par(new=TRUE)
plot(ch.111$expf ~ ch.111$age, 
     type="l", col="darkgray",
     ylim=c(0,350),
     xlim=c(0,100),
     main=NA,
     ylab=NA, 
     xlab=NA)
## par(new=TRUE)
## plot(ch.531$expf ~ ch.531$age, 
##      type="l", col="darkgray",
##      ylim=c(0,350),
##      xlim=c(0,100),
##      main="Trees per Acre",
##      ylab="Trees/Acre > 7.0 Inches", 
##      xlab="Total Age, in years")
################################################################################
## plot the trees per acre, 7.0 inches dbh and larger
## table three (3). 
## reflects pure even-aged second growth stand
par(new=TRUE)
##plot(rcon.111$tpa.7.plus ~ rcon.111$age, 
##plot(rcon.111$expf ~ rcon.111$age, 
plot(res$expf ~ res$age, 
     type="l", col="black",
     ylim=c(0,350),
     xlim=c(0,100),
     main=NA,
     ylab=NA, 
     xlab=NA)
## par(new=TRUE)
## ##plot(rcon.531$tpa.7.plus ~ rcon.531$age, 
## plot(rcon.531$expf ~ rcon.531$age, 
##      type="l", col="black",
##      ylim=c(0,350),
##      xlim=c(0,100),
##      main=NA,
##      ylab=NA, 
##      xlab=NA)
lines(b201.age, b201.tpa, lty=2, col="darkgray")

################################################################################
## plot the basal area for the two models
################################################################################
plot(ch.111$ba ~ ch.111$age, 
     type="n", 
     ylim=c(0,300),
     xlim=c(0,100),
     col="black",
     main="Basal Area",
     ylab="Ft^2/Acre", 
     xlab="Total Age, in years")
grid(lty=1)
par(new=TRUE)
plot(ch.111$ba ~ ch.111$age, 
     type="l", 
     ylim=c(0,300),
     xlim=c(0,100),
     col="darkgray",
     main=NA,
     ylab=NA, 
     xlab=NA)
par(new=TRUE)
## plot(ch.531$ba ~ ch.531$age, 
##      type="l", 
##      ylim=c(0,300),
##      xlim=c(0,100),
##      col="darkgray",
##      main=NA,
##      ylab=NA, 
##      xlab=NA)
## par(new=TRUE)
##plot(rcon.111$ba ~ rcon.111$age, 
plot(res$ba ~ res$age, 
     type="l", 
     ylim=c(0,300),
     xlim=c(0,100),
     col="black",
     main=NA,
     ylab=NA, 
     xlab=NA)
## par(new=TRUE)
## plot(rcon.531$ba ~ rcon.531$age, 
##      type="l", 
##      ylim=c(0,300),
##      xlim=c(0,100),
##      col="black",
##      main=NA,
##      ylab=NA, 
##      xlab=NA)
lines(b201.age, b201.ba, lty=2, col="darkgray")


################################################################################
## plot the quadratic mean diameter for the two models
################################################################################
## plot(ch.111$qmd ~ ch.111$age, 
##      type="n", 
##      ylim=c(0,25),
##      xlim=c(0,100),
##      col="black",
##      main="Mean Stem Diameter",
##      ylab="Inches", 
##      xlab="Total Age, in years")
## grid(lty=1)
## par(new=TRUE)
## plot(ch.111$qmd ~ ch.111$age, 
##      type="l", 
##      ylim=c(0,25),
##      xlim=c(0,100),
##      col="darkgray",
##      main=NA,
##      ylab=NA, 
##      xlab=NA)
## par(new=TRUE)
## plot(ch.531$qmd ~ ch.531$age, 
##      type="l", 
##      ylim=c(0,25),
##      xlim=c(0,100),
##      col="darkgray",
##      main=NA,
##      ylab=NA, 
##      xlab=NA)
## par(new=TRUE)
## plot(rcon.111$qmd ~ rcon.111$age, 
##      type="l", 
##      ylim=c(0,25),
##      xlim=c(0,100),
##      col="black",
##      main=NA,
##      ylab=NA, 
##      xlab=NA)
## par(new=TRUE)
## plot(rcon.531$qmd ~ rcon.531$age, 
##      type="l", 
##      ylim=c(0,25),
##      xlim=c(0,100),
##      col="black",
##      main=NA,
##      ylab=NA, 
##      xlab=NA)



###################################################
### code chunk number 22: simulations.rnw:1099-1114
###################################################
bcmof.diDBH <- function(dbh, tht, cr, hi) {
  p <- 0.25    
  dib <- 1.02453 * dbh^0.88809 * 1.00035^dbh
  X <- (1.0 - sqrt(hi / tht)) / (1.0 - sqrt(p))
  Z = hi / tht
  a = 0.95086 * Z * Z;
  b = -0.18090 * log(Z + 0.001);
  c = 0.61407 * sqrt(Z) +  -0.35105 * exp(Z);
  d = 0.05686 * (dbh / tht); 
  retval <- (dib * X^(a + b + c + d))
  retval[tht < hi] <- 0
  retval
}




###################################################
### code chunk number 23: simulations.rnw:1133-1140
###################################################

merch.height.func <- function(hi, dbh, tht, cr, md) {
  dib <- bcmof.diDBH(dbh, tht, cr, hi)
  diff <- dib - md 
  diff
}



###################################################
### code chunk number 24: simulations.rnw:1150-1160
###################################################

mh <- uniroot(merch.height.func,
                c(0, 27),
                dbh = 45.0,
                tht = 27,
                cr = 0.60,
                md = 10.0)

mh



###################################################
### code chunk number 25: simulations.rnw:1187-1191
###################################################

min.dib.check <- bcmof.diDBH(45, 27, 0.60, mh$root)
min.dib.check



###################################################
### code chunk number 26: simulations.rnw:1227-1243
###################################################

smal.vol <- function(d1, d2, l) {
  c <- 0.0001570796 ## in m^3
  smal.vol <- (0.25 * d1^2 + 0.25 * d2^2) * l * c
  smal.vol
}


log.breaks=c(2,5,12,18,32,999)
log.grades=c("pulp","s4","s3","s2","s1","peeler")
grade.names <- c("Pulp","#4 Sawlog","#3 Sawlog",
                 "#2 Sawlog","#1 Sawlog", "Peeler")

## log.breaks <- c(5,15,30,75,999)
## log.grades <- c("pulp","s3","s2","s1","peeler")



###################################################
### code chunk number 27: log loops
###################################################

mh.bks <- rep(0, length(log.breaks))
for(i in 1:length(log.breaks)) {      
  dbh <- 45.0; tht <- 27; cr <- 0.60          
  if(log.breaks[i] <= bcmof.diDBH(dbh, tht, cr, 0.0)) {    
    mh <- uniroot(merch.height.func,
                  c(0, tht),
                  dbh = dbh,
                  tht = tht,
                  cr = 0.60,
                  md = log.breaks[i])
    mh.bks[i] <- mh$root
  } else {
    mh.bks[i] <- 0.0
  }
}



###################################################
### code chunk number 28: simulations.rnw:1296-1300
###################################################

log.breaks
mh.bks



###################################################
### code chunk number 29: simulations.rnw:1326-1329
###################################################

source("../../scripts/generate_log_volumes.r")



###################################################
### code chunk number 30: fig-log-bucking
###################################################

hi <- 0:tht
dbh <- rep(45.0,length(hi))
tht <- rep(27,length(hi))
cr <- rep(0.6,length(hi))
stem.tpr <- data.frame(cbind(dbh, tht, cr, hi))
stem.tpr$dib <- NA

#for(i in 1:nrow(stem.tpr)){  
#  stem.tpr[i,]$dib <- bcmof.diDBH(
#                                  stem.tpr[i,]$dbh, 
#                                  stem.tpr[i,]$tht, 
#                                  stem.tpr[i,]$cr, 
#                                  stem.tpr[i,]$hi)
#  
#}

## Use vectorized function!

stem.tpr$dib <- with(stem.tpr, bcmof.diDBH(dbh, tht, cr, hi))  

## plot the stem profile
opar <- par(las=1, cex=0.8)
plot(stem.tpr$dib ~ stem.tpr$hi, type="l", 
     ylab=expression(dib),
     xlab=expression(h[i]),
     main="Example Stem with Crosscuts", xlim=c(-1,35))
abline(v=mh.bks, lty=2, lwd=2)
##abline(v=0.3, lty=2, lwd=2)
abline(h=log.breaks, lty=3)
abline(h=0) 
text(x=mh.bks + 4.5,
     y=log.breaks + 2.5,
     labels=paste(log.breaks, " cm @", round(mh.bks, 2), " m"))
par(opar)



###################################################
### code chunk number 31: fig-log-bucking
###################################################

hi <- 0:tht
dbh <- rep(45.0,length(hi))
tht <- rep(27,length(hi))
cr <- rep(0.6,length(hi))
stem.tpr <- data.frame(cbind(dbh, tht, cr, hi))
stem.tpr$dib <- NA

#for(i in 1:nrow(stem.tpr)){  
#  stem.tpr[i,]$dib <- bcmof.diDBH(
#                                  stem.tpr[i,]$dbh, 
#                                  stem.tpr[i,]$tht, 
#                                  stem.tpr[i,]$cr, 
#                                  stem.tpr[i,]$hi)
#  
#}

## Use vectorized function!

stem.tpr$dib <- with(stem.tpr, bcmof.diDBH(dbh, tht, cr, hi))  

## plot the stem profile
opar <- par(las=1, cex=0.8)
plot(stem.tpr$dib ~ stem.tpr$hi, type="l", 
     ylab=expression(dib),
     xlab=expression(h[i]),
     main="Example Stem with Crosscuts", xlim=c(-1,35))
abline(v=mh.bks, lty=2, lwd=2)
##abline(v=0.3, lty=2, lwd=2)
abline(h=log.breaks, lty=3)
abline(h=0) 
text(x=mh.bks + 4.5,
     y=log.breaks + 2.5,
     labels=paste(log.breaks, " cm @", round(mh.bks, 2), " m"))
par(opar)



###################################################
### code chunk number 32: simulations.rnw:1416-1421
###################################################

## the log.breaks are in cm
sp2 <- sp.sums.2(sample.25, log.breaks/2.54, log.grades)
sp2[,c("sm.vol", log.grades)]



###################################################
### code chunk number 33: summary
###################################################

ch <- chambers.1980(ages = 1:200, site = 120, pnba = 1.0)      
ch.stem <- data.frame(plot = 1,
                         sp.code = "DF", 
                         d6 = NA,
                         dbh = ch[50,]$qmd, 
                         tht = ch[50,]$tht,
                         cr = 0.6,
                         n.stems = 1,
                         expf = ch[50,]$expf,
                         crown.width = NA,
                         errors = 0)
ch.plot <- data.frame(plot = 1,
                      elevation = NA,
                      slope = NA,
                      aspect = NA,
                      whc = NA,
                      map = NA,
                      si30 = 120)

ch.50 <- list(plots = ch.plot, plants = ch.stem, age = 50)
class(ch.50) <- "sample.data"

sh2 <- sp.sums.2(ch.50, log.breaks, log.grades)
sh2[,c("qmd","ba","expf","sm.vol",log.grades)]



###################################################
### code chunk number 34: simulations.rnw:1511-1514
###################################################
dia.obs <- rnorm(n = round(ch[50,]$expf)*10, 
                 mean = ch[50,]$qmd,
                 sd = 0.20*ch[50,]$qmd)


###################################################
### code chunk number 35: generate
###################################################

## generate the plant records
plant.obs <- data.frame(plot = 1,
                        sp.code = "DF",
                        d6 = NA,
                        dbh = dia.obs,
                        tht = 4.5 + exp(7.262195456 + 
                              -5.899759104 * 
                              dia.obs^-0.287207389), 
                        cr = 0.40,
                        n.stems = 1,
                        expf = 1/10,
                        crown.width = NA)

plot.obs <- data.frame(plot = 1,
                       elevation = 1000,
                       slope = 0,
                       aspect = 0,
                       whc = NA,
                       map = NA,
                       si30 = 85.0)

ch.50 <- list(plots = plot.obs, 
              plants = plant.obs, 
              age = ch[50,]$age,
              x0 = NA)
class(ch.50)  <- "sample.data"



###################################################
### code chunk number 36: m.o.i.
###################################################

sh2 <- sp.sums.2(ch.50, log.breaks, log.grades)
sh2[, c("qmd", "ba", "expf", "sm.vol", log.grades)]



###################################################
### code chunk number 37: simulations.rnw:1581-1584 (eval = FALSE)
###################################################
## 
## dyn.unload("chambers-1980.so")
## 


###################################################
### code chunk number 38: simulations.rnw:1607-1646
###################################################

## this is where the results get loaded from the earlier version.
load("res")

rcon.res <- res$rconifers
rcon <- rcon.res[rcon.res$veg.ctrl == 1.0, ]
a <- rcon[1:169,]

## with the chambers simulation for the same run
chambers.res <- res$chambers
rch <- chambers.res[chambers.res$veg.ctrl == 1.0,]
b <- rch[1:169,]

## load b.prime here
load(file="b.prime")

b.prime$mai <- b.prime$sm.vol / b.prime$age
## hack to extend b.prime to the starting age.
b.prime <-  rbind(b[1:10,], b.prime[2:nrow(b.prime),])

## these will be a little off. 
## so you can't create a vector with 100 years.
## you have to go from age ? to age ?
alpha <- rep(1, length(a$age))
#a$age < 20

## alpha, now a T x 1 vector of weights
## where the attributes for the two data.frames
## are combined into the 
ds <- 30
de <- 40
alpha[a$age < ds] <- 1
alpha[a$age >= ds & a$age <= de] <- seq(1, 0, length.out=length(ds:de))
alpha[a$age > de] <- 0

c <- alpha * a + (1 - alpha) * b
c.prime <- alpha * a + (1 - alpha) * b.prime




###################################################
### code chunk number 39: fig-mega-plot
###################################################

##opar <- par(mfrow=c(3,2), las=1, cex=0.8)
##opar <- par(mfrow=c(3,2), las=1, cex=0.8, mar=c(4.1,4.1,4.1,2.1) )

## order is bottom, west, north, east
## $mar
## [1] 5.1 4.1 4.1 2.1
opar <- par(mfrow=c(3,2), las=1, cex=0.8, mar=c(4.1,4.1,3.1,2.1) )

################################################################################
## do this for the tpa... 
plot(a$expf ~ a$age, 
     type="n", xlim=c(0,100), ylim=c(0,350), lty=1,
     main="Stems per Acre",
     xlab="Age, in Years",
     ylab="Stems per Acre")
grid(lty=1)
par(new=TRUE)
plot(a$expf ~ a$age, 
     type="l", xlim=c(0,100), ylim=c(0,350), lty=1,
     main=NA,
     xlab=NA,
     ylab=NA)
##lines(c$age, c$tpa.7.plus, lwd=5, col="gray")
lines(c$age, c$expf, lwd=5, col="gray")
par(new=TRUE)
##plot(rcon$tpa.7.plus ~ rcon$age, 
plot(a$expf ~ a$age, 
     type="l", xlim=c(0,100), ylim=c(0,350), lty=1,
     main=NA,
     xlab=NA,
     ylab=NA)
par(new=TRUE)
plot(b$expf ~ b$age, 
     type="l", xlim=c(0,100), ylim=c(0,350), lty=2,
     main=NA,
     xlab=NA,
     ylab=NA)
lines(b201.age, b201.tpa, col="gray", lty=2, lwd=2)

################################################################################
## do this for the basal area... 
##plot(rcon$ba ~ rcon$age, 
plot(a$ba ~ a$age, 
     type="n", xlim=c(0,100), ylim=c(0,300), lty=1,
     main="Stand Basal Area",
     xlab="Age, in Years",
     ylab="Stand Basal Area, in Feet^2")
grid(lty=1)
par(new=TRUE)
plot(a$ba ~ a$age, 
     type="l", xlim=c(0,100), ylim=c(0,300), lty=1,
     main=NA,
     xlab=NA,
     ylab=NA)
par(new=TRUE)
lines(c$age, c$ba, lwd=5, col="gray")
par(new=TRUE)
plot(a$ba ~ a$age, 
     type="l", xlim=c(0,100), ylim=c(0,300), lty=1,
     main=NA,
     xlab=NA,
     ylab=NA)
par(new=TRUE)
plot(b$ba ~ b$age, 
     type="l", xlim=c(0,100), ylim=c(0,300), lty=2,
     main=NA,
     xlab=NA,
     ylab=NA)
lines(b201.age, b201.ba, col="gray", lty=2, lwd=2)
##lines(b210.age, b210.ba, col="gray", lty=2)
##lines(b210.age, b210.tpa, col="gray", lty=2)


################################################################################
## do this for the height... 
plot(a$tht ~ a$age, 
     type="n", xlim=c(0,150), ylim=c(0,250), lty=1,
     main="Mean Stand Height",
     xlab="Age, in Years",
     ylab=expression(feet))
grid(lty=1)
par(new=TRUE)
plot(a$tht ~ a$age, 
     type="l", xlim=c(0,150), ylim=c(0,250), lty=1,
     main=NA,
     xlab=NA,
     ylab=NA)
par(new=TRUE)
lines(c$age, c$tht, lwd=5, col="gray")
par(new=TRUE)
plot(a$tht ~ a$age, 
     type="l", xlim=c(0,150), ylim=c(0,250), lty=1,
     main=NA,
     xlab=NA,
     ylab=NA)
par(new=TRUE)
plot(b$tht ~ b$age, 
     type="l", xlim=c(0,150), ylim=c(0,250), lty=2,
     main=NA,
     xlab=NA,
     ylab=NA)
## peaks at 119 years.
##abline(v=119)
##lines(b201.age, b201.tht, lwd=5, col="gray")
lines(b201.age, b201.tht, col="gray", lty=2, lwd=2)
       

################################################################################
## do this for volume... 
plot(a$sm.vol ~ a$age, 
     type="n", xlim=c(0,150), ylim=c(0,25000), lty=1,
     main="Smalian Volume",
     xlab="Age, in Years",
     ylab=expression(ft^3))
grid(lty=1)
par(new=TRUE)
plot(a$sm.vol ~ a$age, 
     type="l", xlim=c(0,150), ylim=c(0,25000), lty=1,
     main=NA,
     xlab=NA,
     ylab=NA)
par(new=TRUE)
##lines(c$age, c$sm.vol, lwd=5, col="gray")
lines(c.prime$age, c.prime$sm.vol, lwd=5, col="gray")
par(new=TRUE)
plot(a$sm.vol ~ a$age, 
     type="l", xlim=c(0,150), ylim=c(0,25000), lty=1,
     main=NA,
     xlab=NA,
     ylab=NA)
par(new=TRUE)
plot(b$sm.vol ~ b$age, 
     type="l", xlim=c(0,150), ylim=c(0,25000), lty=2,
     main=NA,
     xlab=NA,
     ylab=NA)
lines(b201.age, b201.vol, col="gray", lty=2, lwd=2)
lines(b.prime$age, b.prime$sm.vol, col="black", lty=2, lwd=2)
##lines(c.prime$age, c.prime$sm.vol, lwd=5, col="gray")
##lines(b.prime$age, b.prime$sm.vol, col="black", lty=2, lwd=2)


################################################################################
## mai
################################################################################
## plot the unfixed version
plot(a$mai ~ a$age, 
     type="n", xlim=c(0,100), ylim=c(0,300), lty=1,
     ##main="Mean Annual Increment\n(c=alpha a + (1-\alpha)b)",
     ##main="Mean Annual Increment\nWithout DDIST (c)",
     main="Mean Annual Increment\nfrom Stand Summaries",
     xlab="Age, in Years",
     ylab="Cubic Volume/Age")
grid(lty=1)
par(new=TRUE)
plot(a$mai ~ a$age, 
     type="l", xlim=c(0,100), ylim=c(0,300), lty=1,
     main=NA,
     xlab=NA,
     ylab=NA)
par(new=TRUE)
lines(c$age, c$mai, lwd=5, col="gray")
par(new=TRUE)
plot(a$mai ~ a$age, 
     type="l", xlim=c(0,100), ylim=c(0,300), lty=1,
     main=NA,
     xlab=NA,
     ylab=NA)
par(new=TRUE)
plot(b$mai ~ b$age, 
     type="l", xlim=c(0,100), ylim=c(0,300), lty=2,
     main=NA,
     xlab=NA,
     ylab=NA)
lines(b201.age[1:8], b201.mai, lty=2, col="gray", lwd=2)




################################################################################
## plot the fixed version
plot(a$mai ~ a$age, 
     type="n", xlim=c(0,100), ylim=c(0,300), lty=1,
##     main="Mean Annual Increment\n(c.prime=alpha a + (1-\alpha)b.prime)",
     main="Mean Annual Increment\nfrom Tree Lists",
     xlab="Age, in Years",
     ylab="Cubic Volume/Age")
grid(lty=1)
par(new=TRUE)
plot(a$mai ~ a$age, 
     type="l", xlim=c(0,100), ylim=c(0,300), lty=1,
     main=NA,
     xlab=NA,
     ylab=NA)
par(new=TRUE)
##lines(c$age, c$mai, lwd=5, col="gray")
## c.prime is the combined values, with the 
## ddist correction (flewelling correction)
lines(c.prime$age, c.prime$mai, lwd=5, col="gray")
par(new=TRUE)
plot(a$mai ~ a$age, 
     type="l", xlim=c(0,100), ylim=c(0,300), lty=1,
     main=NA,
     xlab=NA,
     ylab=NA)
lines(b201.mai ~ b201.age[1:8], lty=2, lwd=2, col="gray")
par(new=TRUE)
plot(b.prime$mai ~ b.prime$age, 
     type="l", xlim=c(0,100), ylim=c(0,300), lty=2,
     main=NA,
     xlab=NA,
     ylab=NA)

par(opar)



###################################################
### code chunk number 40: fig-mega-plot
###################################################

##opar <- par(mfrow=c(3,2), las=1, cex=0.8)
##opar <- par(mfrow=c(3,2), las=1, cex=0.8, mar=c(4.1,4.1,4.1,2.1) )

## order is bottom, west, north, east
## $mar
## [1] 5.1 4.1 4.1 2.1
opar <- par(mfrow=c(3,2), las=1, cex=0.8, mar=c(4.1,4.1,3.1,2.1) )

################################################################################
## do this for the tpa... 
plot(a$expf ~ a$age, 
     type="n", xlim=c(0,100), ylim=c(0,350), lty=1,
     main="Stems per Acre",
     xlab="Age, in Years",
     ylab="Stems per Acre")
grid(lty=1)
par(new=TRUE)
plot(a$expf ~ a$age, 
     type="l", xlim=c(0,100), ylim=c(0,350), lty=1,
     main=NA,
     xlab=NA,
     ylab=NA)
##lines(c$age, c$tpa.7.plus, lwd=5, col="gray")
lines(c$age, c$expf, lwd=5, col="gray")
par(new=TRUE)
##plot(rcon$tpa.7.plus ~ rcon$age, 
plot(a$expf ~ a$age, 
     type="l", xlim=c(0,100), ylim=c(0,350), lty=1,
     main=NA,
     xlab=NA,
     ylab=NA)
par(new=TRUE)
plot(b$expf ~ b$age, 
     type="l", xlim=c(0,100), ylim=c(0,350), lty=2,
     main=NA,
     xlab=NA,
     ylab=NA)
lines(b201.age, b201.tpa, col="gray", lty=2, lwd=2)

################################################################################
## do this for the basal area... 
##plot(rcon$ba ~ rcon$age, 
plot(a$ba ~ a$age, 
     type="n", xlim=c(0,100), ylim=c(0,300), lty=1,
     main="Stand Basal Area",
     xlab="Age, in Years",
     ylab="Stand Basal Area, in Feet^2")
grid(lty=1)
par(new=TRUE)
plot(a$ba ~ a$age, 
     type="l", xlim=c(0,100), ylim=c(0,300), lty=1,
     main=NA,
     xlab=NA,
     ylab=NA)
par(new=TRUE)
lines(c$age, c$ba, lwd=5, col="gray")
par(new=TRUE)
plot(a$ba ~ a$age, 
     type="l", xlim=c(0,100), ylim=c(0,300), lty=1,
     main=NA,
     xlab=NA,
     ylab=NA)
par(new=TRUE)
plot(b$ba ~ b$age, 
     type="l", xlim=c(0,100), ylim=c(0,300), lty=2,
     main=NA,
     xlab=NA,
     ylab=NA)
lines(b201.age, b201.ba, col="gray", lty=2, lwd=2)
##lines(b210.age, b210.ba, col="gray", lty=2)
##lines(b210.age, b210.tpa, col="gray", lty=2)


################################################################################
## do this for the height... 
plot(a$tht ~ a$age, 
     type="n", xlim=c(0,150), ylim=c(0,250), lty=1,
     main="Mean Stand Height",
     xlab="Age, in Years",
     ylab=expression(feet))
grid(lty=1)
par(new=TRUE)
plot(a$tht ~ a$age, 
     type="l", xlim=c(0,150), ylim=c(0,250), lty=1,
     main=NA,
     xlab=NA,
     ylab=NA)
par(new=TRUE)
lines(c$age, c$tht, lwd=5, col="gray")
par(new=TRUE)
plot(a$tht ~ a$age, 
     type="l", xlim=c(0,150), ylim=c(0,250), lty=1,
     main=NA,
     xlab=NA,
     ylab=NA)
par(new=TRUE)
plot(b$tht ~ b$age, 
     type="l", xlim=c(0,150), ylim=c(0,250), lty=2,
     main=NA,
     xlab=NA,
     ylab=NA)
## peaks at 119 years.
##abline(v=119)
##lines(b201.age, b201.tht, lwd=5, col="gray")
lines(b201.age, b201.tht, col="gray", lty=2, lwd=2)
       

################################################################################
## do this for volume... 
plot(a$sm.vol ~ a$age, 
     type="n", xlim=c(0,150), ylim=c(0,25000), lty=1,
     main="Smalian Volume",
     xlab="Age, in Years",
     ylab=expression(ft^3))
grid(lty=1)
par(new=TRUE)
plot(a$sm.vol ~ a$age, 
     type="l", xlim=c(0,150), ylim=c(0,25000), lty=1,
     main=NA,
     xlab=NA,
     ylab=NA)
par(new=TRUE)
##lines(c$age, c$sm.vol, lwd=5, col="gray")
lines(c.prime$age, c.prime$sm.vol, lwd=5, col="gray")
par(new=TRUE)
plot(a$sm.vol ~ a$age, 
     type="l", xlim=c(0,150), ylim=c(0,25000), lty=1,
     main=NA,
     xlab=NA,
     ylab=NA)
par(new=TRUE)
plot(b$sm.vol ~ b$age, 
     type="l", xlim=c(0,150), ylim=c(0,25000), lty=2,
     main=NA,
     xlab=NA,
     ylab=NA)
lines(b201.age, b201.vol, col="gray", lty=2, lwd=2)
lines(b.prime$age, b.prime$sm.vol, col="black", lty=2, lwd=2)
##lines(c.prime$age, c.prime$sm.vol, lwd=5, col="gray")
##lines(b.prime$age, b.prime$sm.vol, col="black", lty=2, lwd=2)


################################################################################
## mai
################################################################################
## plot the unfixed version
plot(a$mai ~ a$age, 
     type="n", xlim=c(0,100), ylim=c(0,300), lty=1,
     ##main="Mean Annual Increment\n(c=alpha a + (1-\alpha)b)",
     ##main="Mean Annual Increment\nWithout DDIST (c)",
     main="Mean Annual Increment\nfrom Stand Summaries",
     xlab="Age, in Years",
     ylab="Cubic Volume/Age")
grid(lty=1)
par(new=TRUE)
plot(a$mai ~ a$age, 
     type="l", xlim=c(0,100), ylim=c(0,300), lty=1,
     main=NA,
     xlab=NA,
     ylab=NA)
par(new=TRUE)
lines(c$age, c$mai, lwd=5, col="gray")
par(new=TRUE)
plot(a$mai ~ a$age, 
     type="l", xlim=c(0,100), ylim=c(0,300), lty=1,
     main=NA,
     xlab=NA,
     ylab=NA)
par(new=TRUE)
plot(b$mai ~ b$age, 
     type="l", xlim=c(0,100), ylim=c(0,300), lty=2,
     main=NA,
     xlab=NA,
     ylab=NA)
lines(b201.age[1:8], b201.mai, lty=2, col="gray", lwd=2)




################################################################################
## plot the fixed version
plot(a$mai ~ a$age, 
     type="n", xlim=c(0,100), ylim=c(0,300), lty=1,
##     main="Mean Annual Increment\n(c.prime=alpha a + (1-\alpha)b.prime)",
     main="Mean Annual Increment\nfrom Tree Lists",
     xlab="Age, in Years",
     ylab="Cubic Volume/Age")
grid(lty=1)
par(new=TRUE)
plot(a$mai ~ a$age, 
     type="l", xlim=c(0,100), ylim=c(0,300), lty=1,
     main=NA,
     xlab=NA,
     ylab=NA)
par(new=TRUE)
##lines(c$age, c$mai, lwd=5, col="gray")
## c.prime is the combined values, with the 
## ddist correction (flewelling correction)
lines(c.prime$age, c.prime$mai, lwd=5, col="gray")
par(new=TRUE)
plot(a$mai ~ a$age, 
     type="l", xlim=c(0,100), ylim=c(0,300), lty=1,
     main=NA,
     xlab=NA,
     ylab=NA)
lines(b201.mai ~ b201.age[1:8], lty=2, lwd=2, col="gray")
par(new=TRUE)
plot(b.prime$mai ~ b.prime$age, 
     type="l", xlim=c(0,100), ylim=c(0,300), lty=2,
     main=NA,
     xlab=NA,
     ylab=NA)

par(opar)



###################################################
### code chunk number 41: simulations.rnw:1980-1997
###################################################

a.max.mai <- a[which.max(a$mai),]
b.max.mai <- b[which.max(b$mai),]
b.prime.max.mai <- b.prime[b.prime$mai == max(b.prime$mai),]

c.max.mai <- c[c$mai == max(c$mai),]
c.max.mai$col <- "c"

c.prime.max.mai <- c.prime[c.prime$mai == max(c.prime$mai),]
c.prime.max.mai$col <- "c.prime"

max.mais <- rbind(a.max.mai, 
                  b.max.mai,
                  b.prime.max.mai,
                  c.max.mai,
                  c.prime.max.mai)



###################################################
### code chunk number 42: simulations.rnw:2007-2011
###################################################

small.mai <- max.mais[max.mais$mai == min(max.mais$mai) & !is.na(max.mais$col),]
biggest.mai <- max.mais[max.mais$mai == max(max.mais$mai) & !is.na(max.mais$col),]



###################################################
### code chunk number 43: simulations.rnw:2030-2040
###################################################

cap <- paste("Statistics for simulations at the maximum mean annual increment for each of the models.") 
basic <- c("col", "age","qmd","tht","ba","expf", "tpa.7.plus", "sm.vol")
mm2 <- max.mais[,basic]

names( mm2 ) <- c("Simulation", "Age","QMD (in)","Height (ft)","Basal Area ($\\mbox{ft}^2$)", "TPA", "TPA 7\'\'+", "Volume ($\\mbox{ft}^3$)")

latex(mm2, rowlabel=NULL, rowname=NULL, file="",caption=cap,
     label="tab:max-mais-metrics", digits=1, placement="!tp")



###################################################
### code chunk number 44: fig-merging-two-models-soln-vol-dist
###################################################

## you have to reduce the vertical space between the images.
##opar <- par(mfrow=c(3,2), las=1, cex=0.8, mar=c(0,0,0,0), oma=c(0,0,0,0) )

##opar <- par(mfrow=c(3,2), las=1, cex=0.8, mar=c(.5,.5,.5,.5), oma=c(0,0,0,0) )

## order is bottom, west, north, east
## $mar
## [1] 5.1 4.1 4.1 2.1
##opar <- par(mfrow=c(3,2), las=1, cex=0.8, mar=c(4.1,4.1,4.1,2.1) )
opar <- par(mfrow=c(3,2), las=1, cex=0.8, mar=c(4.1,4.1,3.1,2.1) )


grd <- 1
ylim <- c(0,50)
plot(a[,c(log.grades[grd])] ~ a[,c("age")], main=grade.names[grd], type="n", lty=1, ylim=ylim, xlim=c(0,150), xlab="Age, in Years", ylab=expression(feet^3))
grid(lty=1)
par(new=TRUE)
plot(c.prime[,c(log.grades[grd])] ~ c.prime[,c("age")], main=NA, lwd=5, lty=1, col="gray", type="l", ylim=ylim, xlim=c(0,150), xlab=NA, ylab=NA)
par(new=TRUE)
## a is a thin, solid, black line
plot(a[,c(log.grades[grd])] ~ a[,c("age")], main=NA, lwd=2, lty=1, col="black", type="l", ylim=ylim, xlim=c(0,150), xlab="Age, in Years", ylab=expression(feet^3))
par(new=TRUE)
## b is a thin, solid, gray line
plot(b[,c(log.grades[grd])] ~ b[,c("age")], main=NA, lwd=2, lty=1, col="gray", type="l", ylim=ylim, xlim=c(0,150), xlab=NA, ylab=NA)
par(new=TRUE)
## c' is a thin dashed black line, over the thick solid gray line
plot(c.prime[,c(log.grades[grd])] ~ c.prime[,c("age")], main=NA, lwd=2, lty=2, col="black", type="l", ylim=ylim, xlim=c(0,150), xlab=NA, ylab=NA)

################################################################################
## plot the #2 log grade results
## plot(a[,c(log.grades[2])] ~ a[,c("age")], main=grade.names[2], type="n", lty=1, ylim=c(0,1000), xlim=c(0,100), xlab="Age, in Years", ylab=expression(feet^3))
## grid(lty=1)
## par(new=TRUE)
## plot(c.prime[,c(log.grades[2])] ~ c.prime[,c("age")], main=NA,  lwd=5, lty=1,col="gray", type="l", ylim=c(0,1000), xlim=c(0,100), xlab=NA, ylab=NA)
## par(new=TRUE)
## ## a is a thin solid black line
## plot(a[,c(log.grades[2])] ~ a[,c("age")], main=NA, lwd=2, lty=1, col="black", type="l", ylim=c(0,1000), xlim=c(0,100), xlab="Age, in Years", ylab=expression(feet^3))
## par(new=TRUE)
## ## b is a thin solid gray line
## plot(b[,c(log.grades[2])] ~ b[,c("age")], main=NA, lwd=2, lty=1, col="gray", type="l", ylim=c(0,1000), xlim=c(0,100), xlab=NA, ylab=NA)
## ##par(new=TRUE)
## ##plot(c[,c(log.grades[2])] ~ c[,c("age")], main=NA, type="l", lty=2, col="gray", ylim=c(0,1000), xlim=c(0,100), xlab=NA, ylab=NA)
## par(new=TRUE)
## ## c.prime is a thin dashed black line
## plot(c.prime[,c(log.grades[2])] ~ c.prime[,c("age")], main=NA, lwd=2, lty=2, col="black", type="l", ylim=c(0,1000), xlim=c(0,100), xlab=NA, ylab=NA)

grd <- 2
ylim <- c(0,1000)
plot(a[,c(log.grades[grd])] ~ a[,c("age")], main=grade.names[grd], type="n", lty=1, ylim=ylim, xlim=c(0,150), xlab="Age, in Years", ylab=expression(feet^3))
grid(lty=1)
par(new=TRUE)
plot(c.prime[,c(log.grades[grd])] ~ c.prime[,c("age")], main=NA, lwd=5, lty=1, col="gray", type="l", ylim=ylim, xlim=c(0,150), xlab=NA, ylab=NA)
par(new=TRUE)
## a is a thin, solid, black line
plot(a[,c(log.grades[grd])] ~ a[,c("age")], main=NA, lwd=2, lty=1, col="black", type="l", ylim=ylim, xlim=c(0,150), xlab="Age, in Years", ylab=expression(feet^3))
par(new=TRUE)
## b is a thin, solid, gray line
plot(b[,c(log.grades[grd])] ~ b[,c("age")], main=NA, lwd=2, lty=1, col="gray", type="l", ylim=ylim, xlim=c(0,150), xlab=NA, ylab=NA)
par(new=TRUE)
## c' is a thin dashed black line, over the thick solid gray line
plot(c.prime[,c(log.grades[grd])] ~ c.prime[,c("age")], main=NA, lwd=2, lty=2, col="black", type="l", ylim=ylim, xlim=c(0,150), xlab=NA, ylab=NA)

################################################################################
## plot the #3 log grade results
ylim <- c(0,15000)
plot(a[,c(log.grades[3])] ~ a[,c("age")], main=grade.names[3], type="n", lty=1, ylim=ylim, xlim=c(0,100), xlab="Age, in Years", ylab=expression(feet^3))
grid(lty=1)
par(new=TRUE)
## c.prime is wide, gray, solid
plot(c.prime[,c(log.grades[3])] ~ c.prime[,c("age")], main=NA, type="l",lwd=5, col="gray",  lty=1, ylim=ylim, xlim=c(0,100), xlab=NA, ylab=NA)
par(new=TRUE)
## a is the thin, black, solid
plot(a[,c(log.grades[3])] ~ a[,c("age")], main=NA, type="l", lwd=2, col="black", lty=1, ylim=ylim, xlim=c(0,100), xlab="Age, in Years", ylab=expression(feet^3))
par(new=TRUE)
## b is thin, gray, solid
plot(b[,c(log.grades[3])] ~ b[,c("age")], main=NA, type="l", lwd=2, col="gray",lty=1,  ylim=ylim, xlim=c(0,100), xlab=NA, ylab=NA)
par(new=TRUE)
## c is thin, gray, dashed
#plot(c[,c(log.grades[3])] ~ c[,c("age")], main=NA, type="l", lwd=1, col="gray", lty=2,  ylim=ylim, xlim=c(0,100), xlab=NA, ylab=NA)
#par(new=TRUE)
## c.prime thin, black, dashed
plot(c.prime[,c(log.grades[3])] ~ c.prime[,c("age")], main=NA, type="l",  lwd=2, col="black", lty=2, ylim=ylim, xlim=c(0,100), xlab=NA, ylab=NA)


## a = thin solid black
## b = thick dashed gray
## c = thin dashed gray
## c.prime = thin dashed black
  
################################################################################
## plot the #4 log grade results
grd <- 4
ylim <- c(0,15000)
plot(a[,c(log.grades[grd])] ~ a[,c("age")], main=grade.names[grd], type="n", lty=1, ylim=ylim, xlim=c(0,150), xlab="Age, in Years", ylab=expression(feet^3))
grid(lty=1)
par(new=TRUE)
plot(c.prime[,c(log.grades[grd])] ~ c.prime[,c("age")], main=NA, lwd=5, lty=1, col="gray", type="l", ylim=ylim, xlim=c(0,150), xlab=NA, ylab=NA)
par(new=TRUE)
## a is a thin, solid, black line
plot(a[,c(log.grades[grd])] ~ a[,c("age")], main=NA, lwd=2, lty=1, col="black", type="l", ylim=ylim, xlim=c(0,150), xlab="Age, in Years", ylab=expression(feet^3))
par(new=TRUE)
## b is a thin, solid, gray line
plot(b[,c(log.grades[grd])] ~ b[,c("age")], main=NA, lwd=2, lty=1, col="gray", type="l", ylim=ylim, xlim=c(0,150), xlab=NA, ylab=NA)
par(new=TRUE)
## c' is a thin dashed black line, over the thick solid gray line
plot(c.prime[,c(log.grades[grd])] ~ c.prime[,c("age")], main=NA, lwd=2, lty=2, col="black", type="l", ylim=ylim, xlim=c(0,150), xlab=NA, ylab=NA)

grd <- 5
ylim <- c(0,5000)
plot(a[,c(log.grades[grd])] ~ a[,c("age")], main=grade.names[grd], type="n", lty=1, ylim=ylim, xlim=c(0,150), xlab="Age, in Years", ylab=expression(feet^3))
grid(lty=1)
par(new=TRUE)
plot(c.prime[,c(log.grades[grd])] ~ c.prime[,c("age")], main=NA, lwd=5, lty=1, col="gray", type="l", ylim=ylim, xlim=c(0,150), xlab=NA, ylab=NA)
par(new=TRUE)
## a is a thin, solid, black line
plot(a[,c(log.grades[grd])] ~ a[,c("age")], main=NA, lwd=2, lty=1, col="black", type="l", ylim=ylim, xlim=c(0,150), xlab="Age, in Years", ylab=expression(feet^3))
par(new=TRUE)
## b is a thin, solid, gray line
plot(b[,c(log.grades[grd])] ~ b[,c("age")], main=NA, lwd=2, lty=1, col="gray", type="l", ylim=ylim, xlim=c(0,150), xlab=NA, ylab=NA)
par(new=TRUE)
## c' is a thin dashed black line, over the thick solid gray line
plot(c.prime[,c(log.grades[grd])] ~ c.prime[,c("age")], main=NA, lwd=2, lty=2, col="black", type="l", ylim=ylim, xlim=c(0,150), xlab=NA, ylab=NA)


################################################################################
## plot the #6 log grade results
grd <- 6
ylim <- c(0,500)
plot(a[,c(log.grades[grd])] ~ a[,c("age")], main=grade.names[grd], type="n", lty=1, ylim=ylim, xlim=c(0,150), xlab="Age, in Years", ylab=expression(feet^3))
grid(lty=1)
par(new=TRUE)
plot(c.prime[,c(log.grades[grd])] ~ c.prime[,c("age")], main=NA, lwd=5, lty=1, col="gray", type="l", ylim=ylim, xlim=c(0,150), xlab=NA, ylab=NA)
par(new=TRUE)
## a is a thin, solid, black line
plot(a[,c(log.grades[grd])] ~ a[,c("age")], main=NA, lwd=2, lty=1, col="black", type="l", ylim=ylim, xlim=c(0,150), xlab="Age, in Years", ylab=expression(feet^3))
par(new=TRUE)
## b is a thin, solid, gray line
plot(b[,c(log.grades[grd])] ~ b[,c("age")], main=NA, lwd=2, lty=1, col="gray", type="l", ylim=ylim, xlim=c(0,150), xlab=NA, ylab=NA)
##par(new=TRUE)
par(new=TRUE)
## c' is a thin dashed black line, over the thick solid gray line
plot(c.prime[,c(log.grades[grd])] ~ c.prime[,c("age")], main=NA, lwd=2, lty=2, col="black", type="l", ylim=ylim, xlim=c(0,150), xlab=NA, ylab=NA)

## c the a thin, solid, gray line
##plot(c[,c(log.grades[grd])] ~ c[,c("age")], main=NA, lwd=2, lty=1, col="gray", type="l", ylim=ylim, xlim=c(0,150), xlab=NA, ylab=NA)
##par(new=TRUE)

par(opar)




###################################################
### code chunk number 45: fig-merging-two-models-soln-vol-dist
###################################################

## you have to reduce the vertical space between the images.
##opar <- par(mfrow=c(3,2), las=1, cex=0.8, mar=c(0,0,0,0), oma=c(0,0,0,0) )

##opar <- par(mfrow=c(3,2), las=1, cex=0.8, mar=c(.5,.5,.5,.5), oma=c(0,0,0,0) )

## order is bottom, west, north, east
## $mar
## [1] 5.1 4.1 4.1 2.1
##opar <- par(mfrow=c(3,2), las=1, cex=0.8, mar=c(4.1,4.1,4.1,2.1) )
opar <- par(mfrow=c(3,2), las=1, cex=0.8, mar=c(4.1,4.1,3.1,2.1) )


grd <- 1
ylim <- c(0,50)
plot(a[,c(log.grades[grd])] ~ a[,c("age")], main=grade.names[grd], type="n", lty=1, ylim=ylim, xlim=c(0,150), xlab="Age, in Years", ylab=expression(feet^3))
grid(lty=1)
par(new=TRUE)
plot(c.prime[,c(log.grades[grd])] ~ c.prime[,c("age")], main=NA, lwd=5, lty=1, col="gray", type="l", ylim=ylim, xlim=c(0,150), xlab=NA, ylab=NA)
par(new=TRUE)
## a is a thin, solid, black line
plot(a[,c(log.grades[grd])] ~ a[,c("age")], main=NA, lwd=2, lty=1, col="black", type="l", ylim=ylim, xlim=c(0,150), xlab="Age, in Years", ylab=expression(feet^3))
par(new=TRUE)
## b is a thin, solid, gray line
plot(b[,c(log.grades[grd])] ~ b[,c("age")], main=NA, lwd=2, lty=1, col="gray", type="l", ylim=ylim, xlim=c(0,150), xlab=NA, ylab=NA)
par(new=TRUE)
## c' is a thin dashed black line, over the thick solid gray line
plot(c.prime[,c(log.grades[grd])] ~ c.prime[,c("age")], main=NA, lwd=2, lty=2, col="black", type="l", ylim=ylim, xlim=c(0,150), xlab=NA, ylab=NA)

################################################################################
## plot the #2 log grade results
## plot(a[,c(log.grades[2])] ~ a[,c("age")], main=grade.names[2], type="n", lty=1, ylim=c(0,1000), xlim=c(0,100), xlab="Age, in Years", ylab=expression(feet^3))
## grid(lty=1)
## par(new=TRUE)
## plot(c.prime[,c(log.grades[2])] ~ c.prime[,c("age")], main=NA,  lwd=5, lty=1,col="gray", type="l", ylim=c(0,1000), xlim=c(0,100), xlab=NA, ylab=NA)
## par(new=TRUE)
## ## a is a thin solid black line
## plot(a[,c(log.grades[2])] ~ a[,c("age")], main=NA, lwd=2, lty=1, col="black", type="l", ylim=c(0,1000), xlim=c(0,100), xlab="Age, in Years", ylab=expression(feet^3))
## par(new=TRUE)
## ## b is a thin solid gray line
## plot(b[,c(log.grades[2])] ~ b[,c("age")], main=NA, lwd=2, lty=1, col="gray", type="l", ylim=c(0,1000), xlim=c(0,100), xlab=NA, ylab=NA)
## ##par(new=TRUE)
## ##plot(c[,c(log.grades[2])] ~ c[,c("age")], main=NA, type="l", lty=2, col="gray", ylim=c(0,1000), xlim=c(0,100), xlab=NA, ylab=NA)
## par(new=TRUE)
## ## c.prime is a thin dashed black line
## plot(c.prime[,c(log.grades[2])] ~ c.prime[,c("age")], main=NA, lwd=2, lty=2, col="black", type="l", ylim=c(0,1000), xlim=c(0,100), xlab=NA, ylab=NA)

grd <- 2
ylim <- c(0,1000)
plot(a[,c(log.grades[grd])] ~ a[,c("age")], main=grade.names[grd], type="n", lty=1, ylim=ylim, xlim=c(0,150), xlab="Age, in Years", ylab=expression(feet^3))
grid(lty=1)
par(new=TRUE)
plot(c.prime[,c(log.grades[grd])] ~ c.prime[,c("age")], main=NA, lwd=5, lty=1, col="gray", type="l", ylim=ylim, xlim=c(0,150), xlab=NA, ylab=NA)
par(new=TRUE)
## a is a thin, solid, black line
plot(a[,c(log.grades[grd])] ~ a[,c("age")], main=NA, lwd=2, lty=1, col="black", type="l", ylim=ylim, xlim=c(0,150), xlab="Age, in Years", ylab=expression(feet^3))
par(new=TRUE)
## b is a thin, solid, gray line
plot(b[,c(log.grades[grd])] ~ b[,c("age")], main=NA, lwd=2, lty=1, col="gray", type="l", ylim=ylim, xlim=c(0,150), xlab=NA, ylab=NA)
par(new=TRUE)
## c' is a thin dashed black line, over the thick solid gray line
plot(c.prime[,c(log.grades[grd])] ~ c.prime[,c("age")], main=NA, lwd=2, lty=2, col="black", type="l", ylim=ylim, xlim=c(0,150), xlab=NA, ylab=NA)

################################################################################
## plot the #3 log grade results
ylim <- c(0,15000)
plot(a[,c(log.grades[3])] ~ a[,c("age")], main=grade.names[3], type="n", lty=1, ylim=ylim, xlim=c(0,100), xlab="Age, in Years", ylab=expression(feet^3))
grid(lty=1)
par(new=TRUE)
## c.prime is wide, gray, solid
plot(c.prime[,c(log.grades[3])] ~ c.prime[,c("age")], main=NA, type="l",lwd=5, col="gray",  lty=1, ylim=ylim, xlim=c(0,100), xlab=NA, ylab=NA)
par(new=TRUE)
## a is the thin, black, solid
plot(a[,c(log.grades[3])] ~ a[,c("age")], main=NA, type="l", lwd=2, col="black", lty=1, ylim=ylim, xlim=c(0,100), xlab="Age, in Years", ylab=expression(feet^3))
par(new=TRUE)
## b is thin, gray, solid
plot(b[,c(log.grades[3])] ~ b[,c("age")], main=NA, type="l", lwd=2, col="gray",lty=1,  ylim=ylim, xlim=c(0,100), xlab=NA, ylab=NA)
par(new=TRUE)
## c is thin, gray, dashed
#plot(c[,c(log.grades[3])] ~ c[,c("age")], main=NA, type="l", lwd=1, col="gray", lty=2,  ylim=ylim, xlim=c(0,100), xlab=NA, ylab=NA)
#par(new=TRUE)
## c.prime thin, black, dashed
plot(c.prime[,c(log.grades[3])] ~ c.prime[,c("age")], main=NA, type="l",  lwd=2, col="black", lty=2, ylim=ylim, xlim=c(0,100), xlab=NA, ylab=NA)


## a = thin solid black
## b = thick dashed gray
## c = thin dashed gray
## c.prime = thin dashed black
  
################################################################################
## plot the #4 log grade results
grd <- 4
ylim <- c(0,15000)
plot(a[,c(log.grades[grd])] ~ a[,c("age")], main=grade.names[grd], type="n", lty=1, ylim=ylim, xlim=c(0,150), xlab="Age, in Years", ylab=expression(feet^3))
grid(lty=1)
par(new=TRUE)
plot(c.prime[,c(log.grades[grd])] ~ c.prime[,c("age")], main=NA, lwd=5, lty=1, col="gray", type="l", ylim=ylim, xlim=c(0,150), xlab=NA, ylab=NA)
par(new=TRUE)
## a is a thin, solid, black line
plot(a[,c(log.grades[grd])] ~ a[,c("age")], main=NA, lwd=2, lty=1, col="black", type="l", ylim=ylim, xlim=c(0,150), xlab="Age, in Years", ylab=expression(feet^3))
par(new=TRUE)
## b is a thin, solid, gray line
plot(b[,c(log.grades[grd])] ~ b[,c("age")], main=NA, lwd=2, lty=1, col="gray", type="l", ylim=ylim, xlim=c(0,150), xlab=NA, ylab=NA)
par(new=TRUE)
## c' is a thin dashed black line, over the thick solid gray line
plot(c.prime[,c(log.grades[grd])] ~ c.prime[,c("age")], main=NA, lwd=2, lty=2, col="black", type="l", ylim=ylim, xlim=c(0,150), xlab=NA, ylab=NA)

grd <- 5
ylim <- c(0,5000)
plot(a[,c(log.grades[grd])] ~ a[,c("age")], main=grade.names[grd], type="n", lty=1, ylim=ylim, xlim=c(0,150), xlab="Age, in Years", ylab=expression(feet^3))
grid(lty=1)
par(new=TRUE)
plot(c.prime[,c(log.grades[grd])] ~ c.prime[,c("age")], main=NA, lwd=5, lty=1, col="gray", type="l", ylim=ylim, xlim=c(0,150), xlab=NA, ylab=NA)
par(new=TRUE)
## a is a thin, solid, black line
plot(a[,c(log.grades[grd])] ~ a[,c("age")], main=NA, lwd=2, lty=1, col="black", type="l", ylim=ylim, xlim=c(0,150), xlab="Age, in Years", ylab=expression(feet^3))
par(new=TRUE)
## b is a thin, solid, gray line
plot(b[,c(log.grades[grd])] ~ b[,c("age")], main=NA, lwd=2, lty=1, col="gray", type="l", ylim=ylim, xlim=c(0,150), xlab=NA, ylab=NA)
par(new=TRUE)
## c' is a thin dashed black line, over the thick solid gray line
plot(c.prime[,c(log.grades[grd])] ~ c.prime[,c("age")], main=NA, lwd=2, lty=2, col="black", type="l", ylim=ylim, xlim=c(0,150), xlab=NA, ylab=NA)


################################################################################
## plot the #6 log grade results
grd <- 6
ylim <- c(0,500)
plot(a[,c(log.grades[grd])] ~ a[,c("age")], main=grade.names[grd], type="n", lty=1, ylim=ylim, xlim=c(0,150), xlab="Age, in Years", ylab=expression(feet^3))
grid(lty=1)
par(new=TRUE)
plot(c.prime[,c(log.grades[grd])] ~ c.prime[,c("age")], main=NA, lwd=5, lty=1, col="gray", type="l", ylim=ylim, xlim=c(0,150), xlab=NA, ylab=NA)
par(new=TRUE)
## a is a thin, solid, black line
plot(a[,c(log.grades[grd])] ~ a[,c("age")], main=NA, lwd=2, lty=1, col="black", type="l", ylim=ylim, xlim=c(0,150), xlab="Age, in Years", ylab=expression(feet^3))
par(new=TRUE)
## b is a thin, solid, gray line
plot(b[,c(log.grades[grd])] ~ b[,c("age")], main=NA, lwd=2, lty=1, col="gray", type="l", ylim=ylim, xlim=c(0,150), xlab=NA, ylab=NA)
##par(new=TRUE)
par(new=TRUE)
## c' is a thin dashed black line, over the thick solid gray line
plot(c.prime[,c(log.grades[grd])] ~ c.prime[,c("age")], main=NA, lwd=2, lty=2, col="black", type="l", ylim=ylim, xlim=c(0,150), xlab=NA, ylab=NA)

## c the a thin, solid, gray line
##plot(c[,c(log.grades[grd])] ~ c[,c("age")], main=NA, lwd=2, lty=1, col="gray", type="l", ylim=ylim, xlim=c(0,150), xlab=NA, ylab=NA)
##par(new=TRUE)

par(opar)




###################################################
### code chunk number 46: simulations.rnw:2275-2289
###################################################

cap <- paste("Volume metrics for simulations at maximum mean annual increment.") 
volumes <- c("col", "age", log.grades, "sm.vol")
mm3 <- max.mais[,volumes]

grade.names[2:5] <- paste("\\", grade.names[2:5], sep="")

names(mm3) <-  c("Model","Age",grade.names,"Total")

mm3$MAI <- mm3$Total / mm3$Age

latex(mm3, rowlabel=NULL, rowname=NULL, file="",caption=cap, caption.loc="top",
     label="tab:max-mais-volumes", digits=1, where='!h')



###################################################
### code chunk number 47: simulations.rnw:2434-2453
###################################################

hres <- res$rconifers
names(hres)[3] <- "rx" 

sql.file <- file("smc-yields.sql", "w")  
sql.columns <- names(hres)

for(i in 1:nrow(hres)) {
  sql.command <- 
      sprintf("insert into results (%s) values (%s);", 
              paste(as.vector(sql.columns), collapse=","),
              paste(hres[i,], collapse=",")) 
  cat(sql.command, file = sql.file, sep = "\n")
  
}

close(sql.file)




###################################################
### code chunk number 48: simulations.rnw:2678-2680
###################################################
system("rm -fr package-Ch7")
package.skeleton(name = "package-Ch7")


