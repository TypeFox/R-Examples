rm(list=ls())
load("eusilcP.rdata")
load("test0.rdata")

# original population
index <- eusilcP$weights == 1
eusilcP$ageG <- cut(eusilcP$age, 10)
datO <- eusilcP[index,]
tab1a <- table(datO$region,datO$gender, datO$hsize) 
tab1b <- table(datO$region,datO$gender, datO$ecoStat) 
tab1c <- table(datO$region,datO$gender, datO$ageG) 

# pop after sim-annealing
index2 <- unlist(test0$weights) == 1
datN <- eusilcP[index2,]
tab2a <- table(datN$region,datN$gender, datN$hsize)
tab2b <- table(datN$region,datN$gender, datN$ecoStat) 
tab2c <- table(datN$region,datN$gender, datN$ageG)

# region x sex x hsize
pdf(file="male_region_hsize.pdf", width=12, height=6)
par(mfrow=c(1,2))
mosaicplot(tab1a[,1,], main="males x region x hsize (orig)") 
mosaicplot(tab2a[,1,], main="males x region x hsize (new)")
dev.off()
pdf(file="females_region_hsize.pdf", width=12, height=6)
par(mfrow=c(1,2))
mosaicplot(tab1a[,2,], main="females x region x hsize (orig)") 
mosaicplot(tab2a[,2,], main="females x region x hsize (new)")
dev.off()

# region x sex x ecoStat
pdf(file="males_region_ecostat.pdf", width=12, height=6)
par(mfrow=c(1,2))
mosaicplot(tab1b[,1,], main="males x region x ecoStat (orig)")
mosaicplot(tab2b[,1,], main="males x region x ecoStat (new)")
dev.off()
pdf(file="females_region_ecostat.pdf", width=12, height=6)
par(mfrow=c(1,2))
mosaicplot(tab1b[,2,], main="females x region x ecoStat (orig)")
mosaicplot(tab2b[,2,], main="females x region x ecoStat (new)")
dev.off()

# region x sex x ageGroups
pdf(file="males_region_age.pdf", width=12, height=6)
par(mfrow=c(1,2))
mosaicplot(tab1c[,1,], main="males x region x ageG (orig)") # males X region x ecoStat
mosaicplot(tab2c[,1,], main="males x region x ageG (new)") # males X region x ecoStat
dev.off()
pdf(file="females_region_age.pdf", width=12, height=6)
par(mfrow=c(1,2))
mosaicplot(tab1c[,2,], main="females x region x ageG (orig)") # females X region x ecoStat
mosaicplot(tab2c[,2,], main="females x region x ageG (new)") # females X region x ecoStat
dev.off()