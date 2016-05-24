### R code from vignette source 'paper.rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: paper.rnw:39-44
###################################################
library(osc)
data(exampledata)
str(exampledata)
pop.list <- cca(exampledata[,1:2],s=1)
str(pop.list)


###################################################
### code chunk number 2: paper.rnw:49-56
###################################################
  palette(rainbow(12))
  pdf(file="pics/exdat1.pdf", paper="special", width=8, height=4)
  par(mfrow=c(1,2))
  plot(exampledata$x,exampledata$y,col="darkblue", pch=15, xlab="", ylab="", cex=1.2)
  plot(pop.list$cluster$long,pop.list$cluster$lat,col=pop.list$cluster$cluster_id, pch=15, xlab="", ylab="", cex=1.2)
    dev.off()
  palette("default")


###################################################
### code chunk number 3: paper.rnw:69-76
###################################################
#initiate empty matrix
exampledata.pop <- matrix(0, nrow=max(exampledata$x), ncol=max(exampledata$y))

#restructure data
for(i in 1:NROW(exampledata)){
  exampledata.pop[exampledata$x[i],exampledata$y[i]] <- exampledata$z[i]
}


###################################################
### code chunk number 4: paper.rnw:81-83
###################################################
example.result <- cca(exampledata.pop, s=1)
str(example.result)


###################################################
### code chunk number 5: paper.rnw:88-93
###################################################
  pdf(file="pics/exdat2.pdf", paper="special", width=8, height=4)
  par(mfrow=c(1,2))
  image(x=1:20, y=1:20,exampledata.pop, xlab="", ylab="")
  image(x=1:20, y=1:20,example.result, col=c("white",rep(rainbow(12),2)), xlab="", ylab="")
  dev.off()


###################################################
### code chunk number 6: paper.rnw:112-118
###################################################
# create a raster and set the projection information
raster <- raster(extent(0,5,0,5),nrow=5,ncol=5)
raster[c(1,2,3,5,6,10,17,18,22,23,24)] <- 1
proj4string(raster) <- CRS("+proj=longlat")
# get a feeling for the dimensions
summary(distance(raster)[])


###################################################
### code chunk number 7: paper.rnw:121-126
###################################################
# cluster all cells with value 1
# for various cluster distances of 150, 240 and 111 km
cluster <- cca(raster, cell.class=1, s=1.5e+05, unit="m")
cluster2 <- cca(raster, cell.class=1, s=2.4e+05, unit="m")
cluster3 <- cca(raster, cell.class=1, s=1.11e+05, unit="m")


###################################################
### code chunk number 8: paper.rnw:133-134
###################################################
pixel <- cca(raster, cell.class=1, s=1)


###################################################
### code chunk number 9: paper.rnw:136-137
###################################################
str(pixel)


###################################################
### code chunk number 10: paper.rnw:142-151
###################################################
pdf("pics/raster.pdf", width=9, height=3)
par(mfrow=c(1,3))
raster[cellFromXY(raster,cluster$cluster[,1:2])] <- cluster$cluster[,3]
plot(raster,col=rainbow(11),legend=FALSE, main="s = 150 km")
raster[cellFromXY(raster,cluster2$cluster[,1:2])] <- cluster2$cluster[,3]
plot(raster,col=rainbow(11),legend=FALSE, main="s = 240 km")
raster[cellFromXY(raster,cluster3$cluster[,1:2])] <- cluster3$cluster[,3]
plot(raster,col=rainbow(11),legend=FALSE, main="s = 111 km")
dev.off()


###################################################
### code chunk number 11: paper.rnw:164-166
###################################################
data("landcover")
cities <- cca(landcover, cell.class=1, s=2000, unit="m")


###################################################
### code chunk number 12: paper.rnw:171-174
###################################################
str(cities)
result <- landcover*NA
result[cellFromXY(result, cities$cluster[,1:2])] <- cities$cluster[,3]


###################################################
### code chunk number 13: paper.rnw:179-185
###################################################
pdf("pics/landcover.pdf", width=8, height=6)
par(mfrow=c(1,2))
cols <- c("lightblue", "darkred","yellow","orange","green","darkgreen")
plot(landcover, col=cols)
plot(result, col=rainbow(9))
dev.off()


###################################################
### code chunk number 14: paper.rnw:205-234 (eval = FALSE)
###################################################
## library(osc)
## data("population")
## 
## tl <- rep(NA,100)
## tm <- rep(NA,100)
## 
## for(i in 1:100){
##   print(i)
##   tm[i] <- system.time(mat <- cca(population, s=i, mode=3))[[1]]
## }
## 
## index <- which(population[]>0)
## long <- ceiling(index/nrow(population))
## lat <- index - ((long-1)*nrow(population))
## popl <- data.frame(long, lat, cluster_id=rep(0, length(lat)))
## 
## for(i in 1:100){
##   print(i)
##   tl[i] <- system.time(list <- cca(popl,s=i) )[[1]]
## }
## 
## saveRDS(tm, "tm.rds")
## saveRDS(tm, "tl.rds")
## 
## pdf("compare.pdf", width=5, height=4)
## plot(tm, xlab="cluster distance", ylab="time in s", col="darkblue", log="", type="l")
## lines(tl, col="darkred")
## legend("topleft", legend=c("matrix","list"), col=c("darkblue","darkred"), pch="-")
## dev.off()


