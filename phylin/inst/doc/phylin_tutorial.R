### R code from vignette source 'phylin_tutorial.Snw'
### Encoding: UTF-8

###################################################
### code chunk number 1: install (eval = FALSE)
###################################################
## install.packages('phylin') 


###################################################
### code chunk number 2: load
###################################################
library(phylin) 


###################################################
### code chunk number 3: attach
###################################################
data(d.gen)
data(vipers)
data(grid)


###################################################
### code chunk number 4: dgen
###################################################
example(d.gen) 


###################################################
### code chunk number 5: dgen_plot
###################################################
example(d.gen) 


###################################################
### code chunk number 6: hc_plot (eval = FALSE)
###################################################
## hc <- hclust(as.dist(d.gen))
## plot(hc, hang = -1)


###################################################
### code chunk number 7: grid_plot
###################################################
example(grid) 


###################################################
### code chunk number 8: samples_plot
###################################################
example(vipers) 


###################################################
### code chunk number 9: realdist
###################################################
r.dist <- dist(vipers[,1:2])


###################################################
### code chunk number 10: gv
###################################################
gv <- gen.variogram(r.dist, d.gen)


###################################################
### code chunk number 11: multipleVariogram (eval = FALSE)
###################################################
## # Assuming that the names in the tree do not need processing 
## # to correspond to names in the real distance matrix.
## d.gen.multi <- lapply(trees, cophenetic)


###################################################
### code chunk number 12: GV_plot_no_model
###################################################
plot(gv)


###################################################
### code chunk number 13: checkn
###################################################
gv$n


###################################################
### code chunk number 14: summary
###################################################
summary(gv)


###################################################
### code chunk number 15: model1
###################################################
gv <- gv.model(gv) 


###################################################
### code chunk number 16: GV_plot
###################################################
plot(gv)


###################################################
### code chunk number 17: extraGV_plot
###################################################
gv2 <- gv.model(gv, range=8)
gv.linear <- gv.model(gv, model='linear', range=8)
layout(matrix(1:2, 1, 2))
plot(gv2)
plot(gv.linear)


###################################################
### code chunk number 18: summary
###################################################
summary(gv)


###################################################
### code chunk number 19: lin
###################################################
lin <- as.integer(vipers$lin == 1) 
int.krig <- krig(lin, vipers[,1:2], grid, gv, clamp = TRUE)


###################################################
### code chunk number 20: krigImg
###################################################
grid.image(int.krig, grid, main='Kriging with genetic distances', 
           xlab='Longitude', ylab='Latitude', 
           sclab='Lineage interpolation')
points(vipers[,1:2], pch=lin+1) 


###################################################
### code chunk number 21: krigSD
###################################################
grid.image(int.krig, grid, ic='sd', main='Kriging with genetic distances', 
           xlab='Longitude', ylab='Latitude', 
           sclab='Standard Deviation')


###################################################
### code chunk number 22: krigBin
###################################################
lin.krig <- as.integer(int.krig$Z>0.95) 
grid.image(lin.krig, grid, main='Kriging with genetic distances', 
           xlab='Longitude', ylab='Latitude', 
           sclab='Lineage presence')
points(vipers[,1:2], pch=lin+1) 


###################################################
### code chunk number 23: treesholds
###################################################
#regular sampling
regSampling <- seq(0.01, 0.08, 0.005) 

#node sampling (avoiding the tips)
nodeSampling <- hc$height[hc$height > 0.01 & hc$height < max(hc$height)] 

#single threshold
singleSampling <- 0.06

length(regSampling) 
length(nodeSampling)
length(singleSampling)


###################################################
### code chunk number 24: treeSampling
###################################################
layout(matrix(1:3, 1, 3)) 
plot(hc, hang = -1, labels = FALSE, main= 'Regular sampling') 
abline (h=regSampling, col='red', lty=2) 
plot(hc, hang = -1, labels = FALSE, main='Node sampling') 
abline (h=nodeSampling, col='red', lty=2) 
plot(hc, hang = -1, labels = FALSE, main='Single threshold') 
abline (h=singleSampling, col='red', lty=2) 


###################################################
### code chunk number 25: contactcode
###################################################
contact = rep(0, nrow(grid)) # Sums all probabilities 

for (h in regSampling) { 
    lins <- cutree(hc, h=h) 
    print(paste("height =", h, ":", max(lins), "lineages")) #keep track
    ct = rep(1, nrow(grid)) # Product of individual cluster/lineage map 
    for (i in unique(lins)) { 
        lin <- as.integer(lins == i) 
        krg <- krig(lin, vipers[,1:2], grid, gv, clamp = TRUE, verbose=FALSE) 

        # Product of the complement of the cluster ocurrence probability. 
        ct <- ct * (1 - krg$Z)
    } 
    contact = contact + ct 
} 
# Recycle krg with averaged potential contact zones
krg$Z <- contact / length(regSampling) 


###################################################
### code chunk number 26: contactReg
###################################################
grid.image(krg, grid, main='Potential contact zones / Regular sampling', 
           xlab='Longitude', ylab='Latitude', 
           sclab='Prob. of multiple lineage occurrence')
points(vipers[,1:2], cex=0.5)


###################################################
### code chunk number 27: contactNode
###################################################
contact = rep(0, nrow(grid)) # Sums all probabilities 

for (h in nodeSampling) { 
    lins <- cutree(hc, h=h) 
    print(paste("height =", h, ":", max(lins), "lineages")) #keep track
    ct = rep(1, nrow(grid)) # Product of individual cluster/lineage map 
    for (i in unique(lins)) { 
        lin <- as.integer(lins == i) 
        krg <- krig(lin, vipers[,1:2], grid, gv, clamp = TRUE, verbose=FALSE) 

        # Product of the complement of the cluster ocurrence probability. 
        ct <- ct * (1 - krg$Z) 
    } 
    contact = contact + ct 
} 
# Recycle krg with averaged potential contact zones
krg$Z <- contact / length(nodeSampling) 

grid.image(krg, grid, main='Potential contact zones / Node sampling', 
           xlab='Longitude', ylab='Latitude', 
           sclab='Prob. of multiple lineage occurrence')
points(vipers[,1:2], cex=0.5)


###################################################
### code chunk number 28: contactSingle
###################################################
contact = rep(0, nrow(grid)) # Sums all probabilities 

for (h in singleSampling) { 
    lins <- cutree(hc, h=h) 
    print(paste("height =", h, ":", max(lins), "lineages")) #keep track
    ct = rep(1, nrow(grid)) # Product of individual cluster/lineage map 
    for (i in unique(lins)) { 
        lin <- as.integer(lins == i) 
        krg <- krig(lin, vipers[,1:2], grid, gv, clamp = TRUE, verbose=FALSE) 

        # Product of the complement of the cluster ocurrence probability. 
        ct <- ct * (1 - krg$Z)  
    } 
    contact = contact + ct 
} 
# Recycle krg with averaged potential contact zones
krg$Z <- contact / length(singleSampling) 

grid.image(krg, grid, main='Potential contact zones / single threshold', 
           xlab='Longitude', ylab='Latitude', 
           sclab='Prob. of multiple lineage occurrence')
points(vipers[,1:2], cex=0.5)


###################################################
### code chunk number 29: midpoints
###################################################
mp <- midpoints(vipers[,1:2]) 

d.real <- as.matrix(r.dist) #real distances to matrix

fit <- lm(as.vector(d.gen) ~ as.vector(d.real)) 
resid <- matrix(fit$residuals, nrow(vipers), nrow(vipers)) 
dimnames(resid) <- dimnames(d.gen) 
mp$z <- extract.val(resid, mp[,1:2]) 

int <- idw(mp[,5], mp[,3:4], grid)


###################################################
### code chunk number 30: midpoints
###################################################
grid.image(int, grid, main='IDW interpolation', 
        xlab='Longitude', ylab='Latitude', 
        sclab="Residuals of genetic vs. real distances") 

# plot samples connecting lines 
for (i in 1:nrow(mp)) { 
    pair <- as.character(unlist(mp[i,1:2])) 
    x <- c(vipers[pair[1],1], vipers[pair[2],1]) 
    y <- c(vipers[pair[1],2], vipers[pair[2],2]) 
    lines(x, y, lty=2) 
} 

points(vipers[,1:2], pch=16) # plot samples points in black 
points(mp[,3:4], pch=16, col='gray') # plot midpoints in gray


