### R code from vignette source 'multinomial-travel.Rnw'

###################################################
### code chunk number 1: multinomial-travel.Rnw:13-14 (eval = FALSE)
###################################################
## options(width=80)


###################################################
### code chunk number 2: multinomial-travel.Rnw:20-21 (eval = FALSE)
###################################################
## library(mlogit)


###################################################
### code chunk number 3: multinomial-travel.Rnw:27-28 (eval = FALSE)
###################################################
## data(ModeChoice, package="Ecdat")


###################################################
### code chunk number 4: multinomial-travel.Rnw:34-36 (eval = FALSE)
###################################################
## travel.long <- mlogit.data(ModeChoice, choice="mode", shape="long", alt.levels=
## c("air","train","bus","car"))


###################################################
### code chunk number 5: multinomial-travel.Rnw:42-44 (eval = FALSE)
###################################################
## travel.kat.id <- mlogit(mode ~ invt + gc|hinc, data=travel.long)
## summary(travel.kat.id)


###################################################
### code chunk number 6: multinomial-travel.Rnw:49-50 (eval = FALSE)
###################################################
## library(VGAM)


###################################################
### code chunk number 7: multinomial-travel.Rnw:55-73 (eval = FALSE)
###################################################
## travelmode <- matrix(ModeChoice$mode, byrow = T, ncol = 4)
## colnames(travelmode) <- c("air","train","bus","car")
## travelhinc <- matrix(ModeChoice$hinc, byrow = T, ncol = 4)
## travelhinc <- travelhinc[,1]
## travelinvt <- matrix(ModeChoice$invt, byrow = T, ncol = 4)
## colnames(travelinvt) <- c("invtair","invttrain","invtbus","invtcar")
## travelgc <- matrix(ModeChoice$gc, byrow = T, ncol = 4)
## colnames(travelgc) <- c("gcair","gctrain","gcbus","gccar")
## 
## travelinvt <- sweep(travelinvt[,-1], 1, travelinvt[,1])
## travelgc <- sweep(travelgc[,-1], 1, travelgc[,1])
## 
## 
## Invt <- travelinvt[,1]
## Gc <- travelgc[,1]
## 
## traveldat <- cbind(travelhinc, travelinvt, Invt, travelgc, Gc)
## traveldat <- as.data.frame(traveldat)


###################################################
### code chunk number 8: multinomial-travel.Rnw:78-88 (eval = FALSE)
###################################################
## fit <- vglm(travelmode ~ Invt + Gc + travelhinc,
##             multinomial(parallel = FALSE ~ travelhinc, refLevel = 1),
##             xij = list(Invt ~ invttrain + invtbus + invtcar,
##                        Gc ~ gctrain + gcbus + gccar),
##             form2 = ~ Invt + invttrain + invtbus + invtcar +
##                       Gc + gctrain + gcbus + gccar + travelhinc,
##             data = traveldat, trace = TRUE)
## 
## summary(fit)
## summary(travel.kat.id)


###################################################
### code chunk number 9: multinomial-travel.Rnw:93-95 (eval = FALSE)
###################################################
## summary(travel.kat.id)$CoefTable
## summary(fit)@coef3


###################################################
### code chunk number 10: multinomial-travel.Rnw:98-100 (eval = FALSE)
###################################################
## detach(package:VGAM)
## detach(package:mlogit)


