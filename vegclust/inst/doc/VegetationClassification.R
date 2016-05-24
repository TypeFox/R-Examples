### R code from vignette source 'VegetationClassification.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: VegetationClassification.Rnw:20-21
###################################################
options(width=67)


###################################################
### code chunk number 2: VegetationClassification.Rnw:26-27
###################################################
library(vegclust)


###################################################
### code chunk number 3: VegetationClassification.Rnw:32-34
###################################################
data(wetland)
dim(wetland)


###################################################
### code chunk number 4: VegetationClassification.Rnw:37-38
###################################################
wetlandchord = decostand(wetland,"normalize")


###################################################
### code chunk number 5: VegetationClassification.Rnw:41-42
###################################################
dchord = dist(wetlandchord)


###################################################
### code chunk number 6: VegetationClassification.Rnw:177-179
###################################################
wetland.km = vegclust(x = wetlandchord, mobileCenters=3, 
                      method="KM", nstart=20)


###################################################
### code chunk number 7: VegetationClassification.Rnw:182-183
###################################################
names(wetland.km)


###################################################
### code chunk number 8: VegetationClassification.Rnw:186-187
###################################################
t(wetland.km$memb)


###################################################
### code chunk number 9: VegetationClassification.Rnw:190-191
###################################################
round(wetland.km$mobileCenters, dig=3)


###################################################
### code chunk number 10: VegetationClassification.Rnw:194-197
###################################################
wetland.kmdist = vegclustdist(x = dchord, mobileMemb=3, 
                              method="KM", nstart = 20)
names(wetland.kmdist)


###################################################
### code chunk number 11: VegetationClassification.Rnw:200-201
###################################################
wetland.kmdist$mobileCenters


###################################################
### code chunk number 12: VegetationClassification.Rnw:204-205
###################################################
t(wetland.kmdist$memb)


###################################################
### code chunk number 13: VegetationClassification.Rnw:208-210
###################################################
wetland.km$mode
wetland.kmdist$mode


###################################################
### code chunk number 14: VegetationClassification.Rnw:216-217
###################################################
round(t(wetland.km$dist2clusters), dig=2)


###################################################
### code chunk number 15: VegetationClassification.Rnw:220-223
###################################################
wetland.fcm = vegclust(x = wetlandchord, mobileCenters=3, 
                       method="FCM", m=1.2, nstart=20)
round(t(wetland.fcm$memb), dig=3)


###################################################
### code chunk number 16: VegetationClassification.Rnw:228-231
###################################################
groups = defuzzify(wetland.fcm)$cluster
groups
table(groups)


###################################################
### code chunk number 17: VegetationClassification.Rnw:234-237
###################################################
groups = defuzzify(wetland.fcm, method = "cut", alpha = 0.8)$cluster
groups
table(groups, useNA = "always")


###################################################
### code chunk number 18: VegetationClassification.Rnw:240-243
###################################################
wetland.fcm2 = vegclust(x = wetlandchord, mobileCenters=3, 
                       method="FCM", m=10, nstart=20)
round(t(wetland.fcm2$memb), dig=3)


###################################################
### code chunk number 19: VegetationClassification.Rnw:246-248
###################################################
groups2 = defuzzify(wetland.fcm2, method = "cut", alpha = 0.8)$cluster
table(groups2, useNA = "always")


###################################################
### code chunk number 20: VegetationClassification.Rnw:254-257
###################################################
wetland.nc = vegclust(x = wetlandchord, mobileCenters=3,
                       method="NC", m=1.2, dnoise=0.8, nstart=20)
round(t(wetland.nc$memb), dig=2)


###################################################
### code chunk number 21: VegetationClassification.Rnw:260-263
###################################################
groups = defuzzify(wetland.nc)$cluster
groups
table(groups)


###################################################
### code chunk number 22: VegetationClassification.Rnw:266-269
###################################################
groups = defuzzify(wetland.nc, method="cut", alpha=0.8)$cluster
groups
table(groups, useNA = "always")


###################################################
### code chunk number 23: VegetationClassification.Rnw:274-277
###################################################
dist(wetland.km$mobileCenters)
dist(wetland.fcm$mobileCenters)
dist(wetland.nc$mobileCenters)


###################################################
### code chunk number 24: VegetationClassification.Rnw:283-286
###################################################
wetland.kmdd = vegclust(x = wetlandchord, mobileCenters=3, 
                      method="KMdd", nstart=20)
t(wetland.kmdd$memb)


###################################################
### code chunk number 25: VegetationClassification.Rnw:289-290
###################################################
round(wetland.kmdd$mobileCenters, dig=3)


###################################################
### code chunk number 26: VegetationClassification.Rnw:293-296
###################################################
wetland.kmdd = vegclustdist(x = dchord, mobileMemb=3, 
                      method="KMdd", nstart=20)
wetland.kmdd$mobileCenters


###################################################
### code chunk number 27: VegetationClassification.Rnw:301-307
###################################################
wetland.31 = wetlandchord[1:31,]
wetland.31 = wetland.31[,colSums(wetland.31)>0]
dim(wetland.31)
wetland.10 = wetlandchord[-(1:31),]
wetland.10 = wetland.10[,colSums(wetland.10)>0] 
dim(wetland.10)


###################################################
### code chunk number 28: VegetationClassification.Rnw:310-312
###################################################
groups = kmeans(wetland.31, 2)$cluster
groups


###################################################
### code chunk number 29: VegetationClassification.Rnw:315-316
###################################################
wetland.31.km = as.vegclust(wetland.31, groups)


###################################################
### code chunk number 30: VegetationClassification.Rnw:319-320
###################################################
wetland.31.km$method


###################################################
### code chunk number 31: VegetationClassification.Rnw:323-325
###################################################
wetland.10.km = vegclass(wetland.31.km, wetland.10)
defuzzify(wetland.10.km)$cluster


###################################################
### code chunk number 32: VegetationClassification.Rnw:328-332
###################################################
wetland.31.nc = as.vegclust(wetland.31, groups, method="HNC", 
                            dnoise = 0.8)
wetland.10.nc = vegclass(wetland.31.nc, wetland.10)
defuzzify(wetland.10.nc)$cluster


###################################################
### code chunk number 33: VegetationClassification.Rnw:341-346
###################################################
cf = conformveg(wetland.31, wetland.10)
wetland.31.cf<- cf$x
wetland.10.cf<- cf$y
dim(wetland.31.cf)
dim(wetland.10.cf)


###################################################
### code chunk number 34: VegetationClassification.Rnw:352-353
###################################################
fixed = clustcentroid(wetland.31.cf, groups)


###################################################
### code chunk number 35: VegetationClassification.Rnw:359-364
###################################################
wetland.nc = vegclust(wetland.10.cf, mobileCenters=1, 
                      fixedCenters = fixed, 
                      method = wetland.31.nc$method,
                      dnoise=wetland.31.nc$dnoise, nstart=10)
defuzzify(wetland.nc)$cluster


###################################################
### code chunk number 36: VegetationClassification.Rnw:369-374
###################################################
wetland.km = vegclust(wetland.10.cf, mobileCenters=1, 
                      fixedCenters = fixed, 
                      method = "KM",
                      nstart=10)
defuzzify(wetland.km)$cluster


###################################################
### code chunk number 37: VegetationClassification.Rnw:381-386
###################################################
wetland.nc = vegclust(rbind(wetland.31.cf,wetland.10.cf), mobileCenters=1, 
                      fixedCenters = fixed, 
                      method = wetland.31.nc$method,
                      dnoise=wetland.31.nc$dnoise, nstart=10)
defuzzify(wetland.nc)$cluster


###################################################
### code chunk number 38: VegetationClassification.Rnw:394-395
###################################################
groups = c(rep(1, 17), rep(2, 14), rep(3,10))


###################################################
### code chunk number 39: VegetationClassification.Rnw:402-404
###################################################
centroids = clustcentroid(wetlandchord, groups)
round(centroids, dig=3)


###################################################
### code chunk number 40: VegetationClassification.Rnw:408-410
###################################################
medoids = clustmedoid(wetlandchord, groups)
medoids


###################################################
### code chunk number 41: VegetationClassification.Rnw:421-422
###################################################
clustvar(wetlandchord, groups)


###################################################
### code chunk number 42: VegetationClassification.Rnw:434-435
###################################################
clustvar(dchord, groups)


###################################################
### code chunk number 43: VegetationClassification.Rnw:438-439
###################################################
clustvar(wetlandchord)


###################################################
### code chunk number 44: VegetationClassification.Rnw:444-445
###################################################
as.dist(as.matrix(dchord)[medoids,medoids])


###################################################
### code chunk number 45: VegetationClassification.Rnw:449-450
###################################################
dist(centroids)


###################################################
### code chunk number 46: VegetationClassification.Rnw:457-458
###################################################
interclustdist(dchord,groups)


###################################################
### code chunk number 47: VegetationClassification.Rnw:465-466
###################################################
c = clustconst(wetlandchord, memb = as.memb(groups))


###################################################
### code chunk number 48: VegetationClassification.Rnw:470-471
###################################################
d=summary(c, mode="all")


###################################################
### code chunk number 49: VegetationClassification.Rnw:475-476
###################################################
summary(c, mode="cluster", name=names(c)[1])


