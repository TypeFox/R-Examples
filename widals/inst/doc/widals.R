### R code from vignette source 'widals.Snw'

###################################################
### code chunk number 1: widals.Snw:153-156
###################################################
options(width=49)
options(prompt=" ")
options(continue="   ")


###################################################
### code chunk number 2: b1
###################################################
options(stringsAsFactors=FALSE)

k.cpus <- 2 #### set the number of cpus for snowfall

library(widals)
data(O3)

Z.all <- as.matrix(O3$Z)[366:730, ]
locs.all <- O3$locs[ , c(2,1)]
hsa.all <- O3$helevs

xdate <- rownames(Z.all)

tau <- nrow(Z.all)
n.all <- ncol(Z.all)

xgeodesic <- TRUE


###################################################
### code chunk number 3: b2
###################################################
Z <- Z.all
locs <- locs.all
n <- n.all

dateDate <- strptime(xdate, "%Y%m%d")
doy <- as.integer(format(dateDate, "%j"))

Ht <- cbind( sin(2*pi*doy/365), cos(2*pi*doy/365) )

Hs.all <- matrix(1, nrow=n.all)
Hst.ls.all <- NULL


Hs <- Hs.all
Hst.ls <- Hst.ls.all

Ht.original <- Ht

##########################
rm.ndx <- create.rm.ndx.ls(n, 14)
b.lag <- 0

train.rng <- 30:tau 
test.rng <- train.rng

GP <- c(1/10, 1)

k.glob <- 9
rho.upper.limit <- 100
rgr.lower.limit <- 10^(-7)

sds.mx <- seq(2, 0.01, length=k.glob) * matrix(1, k.glob, length(GP))

run.parallel <- TRUE
sfInit(TRUE, k.cpus)

FUN.source <- fun.load.hals.a

FUN.source()

MSS.snow(FUN.source, NA, p.ndx.ls, f.d, sds.mx=sds.mx, 
k.glob, k.loc.coef=7, X = NULL)
sfStop()



###################################################
### code chunk number 4: a2
###################################################

Z.als <- Hals.snow(1, Z = Z, Hs = Hs, Ht = Ht, Hst.ls = Hst.ls, 
b.lag = b.lag, GP.mx = matrix(GP, 1, 2))
resids <- Z-Z.als
sqrt( mean( resids[test.rng, ]^2 ) )



###################################################
### code chunk number 5: a3
###################################################

Hs.all <- cbind(matrix(1, nrow=n.all), hsa.all)
Hs <- Hs.all

GP <- c(1/10, 1)

sfInit(TRUE, k.cpus)
MSS.snow(FUN.source, NA, p.ndx.ls, f.d, sds.mx=sds.mx, 
k.glob, k.loc.coef=7, X = NULL)
sfStop()



###################################################
### code chunk number 6: a4
###################################################

Hst.ls.all <- H.Earth.solar(locs[ , 2], locs[ , 1], dateDate)
Hst.ls <- Hst.ls.all

GP <- c(1/10, 1)

sfInit(TRUE, k.cpus)
MSS.snow(FUN.source, NA, p.ndx.ls, f.d, sds.mx=sds.mx, 
k.glob, k.loc.coef=7, X = NULL)
sfStop()



###################################################
### code chunk number 7: a5
###################################################

Hst.ls.all2 <- list()
for(tt in 1:tau) {
    Hst.ls.all2[[tt]] <- cbind(Hst.ls.all[[tt]], Hst.ls.all[[tt]]*hsa.all)
    colnames(Hst.ls.all2[[tt]]) <- c("ISA", "ISAxElev")
}
Hst.ls <- Hst.ls.all2

GP <- c(1/10, 1)

sfInit(TRUE, k.cpus)
MSS.snow(FUN.source, NA, p.ndx.ls, f.d, sds.mx=sds.mx, 
k.glob, k.loc.coef=7, X = NULL)
sfStop()



###################################################
### code chunk number 8: a6
###################################################

Hals.ses(Z, Hs, Ht, Hst.ls, GP[1], GP[2], b.lag, test.rng)



###################################################
### code chunk number 9: a7
###################################################
Z <- Z.all

Hst.ls <- stnd.Hst.ls(Hst.ls.all2)$sHst.ls
Hs <- stnd.Hs(Hs.all)$sHs
Ht <- stnd.Ht(Ht, nrow(Hs))

GP <- c(1/10, 1)

sfInit(TRUE, k.cpus)
MSS.snow(FUN.source, NA, p.ndx.ls, f.d, sds.mx=sds.mx, 
k.glob, k.loc.coef=7, X = NULL)
sfStop()



###################################################
### code chunk number 10: a8
###################################################

z.mean <- mean(Z.all)
z.sd <- sd(as.vector(Z.all))

Z <- (Z.all - z.mean) / z.sd
GP <- c(1/10, 1)

sfInit(TRUE, k.cpus)
MSS.snow(FUN.source, NA, p.ndx.ls, f.d, sds.mx=sds.mx, 
k.glob, k.loc.coef=7, X = NULL)
sfStop()

our.cost * z.sd



###################################################
### code chunk number 11: a9 (eval = FALSE)
###################################################
## 
## Z <- Z.all
## Hst.ls <- Hst.ls.all2
## Hs <- Hs.all
## Ht <- Ht.original
## 
## FUN.source <- fun.load.widals.a
## 
## d.alpha.lower.limit <- 0
## 
## GP <- c(1/1000, 1, 0.01, 3, 1)
## cv <- 2
## lags <- c(0)
## b.lag <- 0
## 
## sds.mx <- seq(2, 0.01, length=k.glob) * matrix(1, k.glob, length(GP))
## ltco <- -10
## stnd.d <- TRUE
## 
## sfInit(TRUE, k.cpus)
## FUN.source()
## 
## MSS.snow(FUN.source, NA, p.ndx.ls, f.d, sds.mx=sds.mx, 
## k.glob, k.loc.coef=7, X = NULL)
## sfStop()
## 


###################################################
### code chunk number 12: a10
###################################################

rm.ndx <- 1:n

Z <- Z.all
Hst.ls <- Hst.ls.all2
Hs <- Hs.all
Ht <- Ht.original

FUN.source <- fun.load.widals.a

d.alpha.lower.limit <- 0

GP <- c(1/1000, 1, 0.01, 3, 1)
cv <- -2
lags <- c(0)
b.lag <- -1

sds.mx <- seq(2, 0.01, length=k.glob) * matrix(1, k.glob, length(GP))
ltco <- -10
stnd.d <- TRUE

sfInit(TRUE, k.cpus)
FUN.source()

MSS.snow(FUN.source, NA, p.ndx.ls, f.d, sds.mx=sds.mx, 
k.glob, k.loc.coef=7, X = NULL)
sfStop()



###################################################
### code chunk number 13: a10
###################################################

Hals.ses(Z, Hs, Ht, Hst.ls, GP[1], GP[2], b.lag, test.rng)



###################################################
### code chunk number 14: a11
###################################################

Z <- Z.all
Hst.ls <- Hst.ls.all2
Hs <- Hs.all
Ht <- Ht.original


Hs0 <- cbind(matrix(1, length(O3$helevs0)), O3$helevs0)
Hst0isa.ls <- H.Earth.solar(O3$locs0[ , 1], O3$locs0[ , 2], dateDate)
Hst0.ls <- list()
for(tt in 1:tau) {
	Hst0.ls[[tt]] <- cbind(Hst0isa.ls[[tt]], Hst0isa.ls[[tt]]*O3$helevs0)
	colnames(Hst0.ls[[tt]]) <- c("ISA", "ISAxElev")
}

locs0 <- O3$locs0[ , c(2,1)]

Z0.hat <- widals.predict(Z, Hs, Ht, Hst.ls, locs, lags, b.lag, Hs0, Hst0.ls, 
locs0=locs0, geodesic=xgeodesic, wrap.around=NULL, GP, stnd.d, ltco)[10:tau, ]

ydate <- xdate[10:tau]

#xcol.vec <- heat.colors(max(round(Z0.hat)))
#xcol.vec <- rev(rainbow(max(round(Z0.hat))))
xcol.vec <- rev(rainbow(630)[1:max(round(Z0.hat))])

xleg.vals <- round( seq(1, max(Z0.hat)-1, length=(5)) / 1 ) * 1
xleg.cols <- xcol.vec[xleg.vals+1]



###################################################
### code chunk number 15: a11 (eval = FALSE)
###################################################
## for(tt in 1:nrow(Z0.hat)) {
##     plot(0, 0, xlim=c(-124.1, -113.9), ylim=c(32.5, 42), type="n", main=ydate[tt])
##     ## points(locs0[ , c(2,1)], cex=Z0.hat[tt,] / 30) ## uncomment to see sites
##     this.zvec <- round(Z0.hat[tt,])
##     this.zvec[ this.zvec < 1] <- 1
##     this.color <- xcol.vec[ this.zvec ]
##     points(locs0[ , c(2,1)], cex=1.14,  col=this.color, pch=19 )
##     #points(locs[ , c(2,1)], cex=Z[tt,] / 30, col="red")
##     legend(-116, 40, legend=rev(xleg.vals), fill=FALSE, col=rev(xleg.cols), 
##     border=NA, bty="n", text.col=rev(xleg.cols))
## 	Sys.sleep(0.1)
## }


###################################################
### code chunk number 16: widals.Snw:511-523
###################################################
ydate <- xdate[10:tau]
tt <- 180
plot(0, 0, xlim=c(-124.1, -113.9), ylim=c(32.5, 42), type="n", main=ydate[tt], 
xlab="", ylab="")
## points(locs0[ , c(2,1)], cex=Z0.hat[tt,] / 30) ## uncomment to see sites
this.zvec <- round(Z0.hat[tt,])
this.zvec[ this.zvec < 1] <- 1
this.color <- xcol.vec[ this.zvec ]
points(locs0[ , c(2,1)], cex=1.14,  col=this.color, pch=19 )
#points(locs[ , c(2,1)], cex=Z[tt,] / 30, col="red")
legend(-116, 40, legend=rev(xleg.vals), fill=FALSE, col=rev(xleg.cols), border=NA, 
bty="n", text.col=rev(xleg.cols))


