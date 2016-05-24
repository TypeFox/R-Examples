### R code from vignette source 'st.Rnw'

###################################################
### code chunk number 1: st.Rnw:89-93
###################################################
library(spacetime)
rm(list = ls())
data(air)
ls()


###################################################
### code chunk number 2: st.Rnw:102-108
###################################################
if (!exists("rural"))
	rural = STFDF(stations, dates, data.frame(PM10 = as.vector(air)))
rr = rural[,"2005::2010"]
unsel = which(apply(as(rr, "xts"), 2, function(x) all(is.na(x))))
r5to10 = rr[-unsel,]
summary(r5to10)


###################################################
### code chunk number 3: st.Rnw:112-114
###################################################
rn = row.names(r5to10@sp)[4:7]
rn


###################################################
### code chunk number 4: st.Rnw:119-124 (eval = FALSE)
###################################################
## par(mfrow=c(2,2))
## # select 4, 5, 6, 7
## for(i in rn) 
## 	acf(na.omit(r5to10[i,]), main = i)
## par(mfrow=c(1,1))


###################################################
### code chunk number 5: st.Rnw:129-134
###################################################
par(mfrow=c(2,2))
# select 4, 5, 6, 7
rn = row.names(r5to10@sp)[4:7]
for(i in rn) 
	acf(na.omit(r5to10[i,]), main = i)


###################################################
### code chunk number 6: st.Rnw:143-144 (eval = FALSE)
###################################################
## acf(na.omit(as(r5to10[rn,], "xts")))


###################################################
### code chunk number 7: st.Rnw:149-150
###################################################
acf(na.omit(as(r5to10[rn,], "xts")))


###################################################
### code chunk number 8: st.Rnw:174-175 (eval = FALSE)
###################################################
## acf(na.omit(as(r5to10[4:10,], "xts")))


###################################################
### code chunk number 9: st.Rnw:180-182
###################################################
library(sp)
print(spDists(r5to10[4:10,]@sp), digits=3)


###################################################
### code chunk number 10: st.Rnw:189-190
###################################################
rs = sample(dim(r5to10)[2], 100)


###################################################
### code chunk number 11: st.Rnw:195-197
###################################################
lst = lapply(rs, function(i) { x = r5to10[,i]; x$ti = i; rownames(x@coords) = NULL; x} )
pts = do.call(rbind, lst)


###################################################
### code chunk number 12: st.Rnw:200-202
###################################################
library(gstat)
v = variogram(PM10~ti, pts[!is.na(pts$PM10),], dX=0)


###################################################
### code chunk number 13: st.Rnw:205-208 (eval = FALSE)
###################################################
## # plot(v, fit.variogram(v, vgm(1, "Exp", 200, 1)))
## vmod = fit.variogram(v, vgm(100, "Exp", 200))
## plot(v, vmod)


###################################################
### code chunk number 14: st.Rnw:212-215
###################################################
# plot(v, fit.variogram(v, vgm(1, "Exp", 200, 1)))
vmod = fit.variogram(v, vgm(100, "Exp", 200))
print(plot(v, vmod))


###################################################
### code chunk number 15: st.Rnw:223-224
###################################################
vmod


###################################################
### code chunk number 16: st.Rnw:229-230
###################################################
dim(r5to10)


###################################################
### code chunk number 17: st.Rnw:235-236 (eval = FALSE)
###################################################
## vv = variogram(PM10~1, r5to10, width=20, cutoff = 200, tlags=0:5)


###################################################
### code chunk number 18: st.Rnw:241-242 (eval = FALSE)
###################################################
## vv = variogram(PM10~1, r5to10, width=20, cutoff = 200, tlags=0:5)


###################################################
### code chunk number 19: st.Rnw:248-249
###################################################
data(vv)


###################################################
### code chunk number 20: st.Rnw:252-253
###################################################
vv <- vv[c("np", "dist", "gamma", "id", "timelag", "spacelag")]


###################################################
### code chunk number 21: st.Rnw:258-260 (eval = FALSE)
###################################################
## plot(vv)
## plot(vv, map = FALSE)


###################################################
### code chunk number 22: st.Rnw:264-266
###################################################
print(plot(vv), split = c(1,1,1,2), more = TRUE)
print(plot(vv, map = FALSE), split = c(1,2,1,2))


###################################################
### code chunk number 23: st.Rnw:278-283
###################################################
metricVgm <- vgmST("metric",
                   joint=vgm(50,"Exp",100,0),
                   stAni=50)

metricVgm <- fit.StVariogram(vv, metricVgm)


###################################################
### code chunk number 24: st.Rnw:288-289
###################################################
attr(metricVgm, "optim")$value


###################################################
### code chunk number 25: st.Rnw:294-295 (eval = FALSE)
###################################################
## plot(vv, metricVgm)


###################################################
### code chunk number 26: st.Rnw:300-301
###################################################
print(plot(vv, metricVgm))


###################################################
### code chunk number 27: st.Rnw:311-319
###################################################
sepVgm <- vgmST("separable",
                space=vgm(0.9,"Exp", 123, 0.1),
                time =vgm(0.9,"Exp", 2.9, 0.1),
                sill=100)

sepVgm <- fit.StVariogram(vv, sepVgm, method = "L-BFGS-B",
                          lower = c(10,0,0.01,0,1),
                          upper = c(500,1,20,1,200))


###################################################
### code chunk number 28: st.Rnw:325-327
###################################################
attr(sepVgm, "optim")$value
plot(vv, list(sepVgm, metricVgm))


###################################################
### code chunk number 29: st.Rnw:332-333
###################################################
print(plot(vv, list(sepVgm, metricVgm)))


###################################################
### code chunk number 30: st.Rnw:341-347 (eval = FALSE)
###################################################
## library(lattice)
## plot(vv, list(sepVgm, metricVgm), all=T, wireframe=T, zlim=c(0,120),
##      zlab=NULL,
##      xlab=list("distance (km)", rot=30),
##      ylab=list("time lag (days)", rot=-35),
##      scales=list(arrows=F, z = list(distance = 5)))


###################################################
### code chunk number 31: st.Rnw:358-364
###################################################
library(lattice)
print(plot(vv, list(sepVgm, metricVgm), all=T, wireframe=T, zlim=c(0,120),
           zlab=NULL,
           xlab=list("distance (km)", rot=30),
           ylab=list("time lag (days)", rot=-35),
           scales=list(arrows=F, z = list(distance = 5))))


###################################################
### code chunk number 32: st.Rnw:385-386 (eval = FALSE)
###################################################
## demo(gstat3D)


