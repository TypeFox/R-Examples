PharmPow_crossover <- function(datafullA,
                          datareducedA,
                          datafullB,
                          datareducedB,
                          maxIDsA=200,
                          maxIDsB=200,
                          nresamples=10000,
                          power=80,
                          designAlab="Number of patients sampling design A",
                          designBlab="Number of patients sampling design B"){

### read in data

datafullA      <- read.table(datafullA,header=TRUE,skip=1)
datareducedA   <- read.table(datareducedA,header=TRUE,skip=1)
datafullB      <- read.table(datafullB,header=TRUE,skip=1)
datareducedB   <- read.table(datareducedB,header=TRUE,skip=1)

Ninda   <- seq(0,maxIDsA,by=1)
Nindb   <- seq(0,maxIDsB,by=1)

Nsam    <- seq(1,nresamples,by=1)

### construct sample dataset###

dvector1   <- c(datafullA$ID,datafullA$OBJ,datareducedA$OBJ)
dA1        <- matrix(dvector1,length(datafullA$ID),3)
dA1        <- as.data.frame(dA1)
names(dA1) <- c("ID","OFV_FULL","OFV_RED")
dA1$FLAG   <- 0

dvector2   <- c(datafullB$ID,datafullB$OBJ,datareducedB$OBJ)
dA2        <- matrix(dvector2,length(datafullB$ID),3)
dA2        <- as.data.frame(dA2)
names(dA2) <- c("ID","OFV_FULL","OFV_RED")
dA2$FLAG   <- 1

d        <- rbind(dA1,dA2)

d$iOFV   <- d$OFV_FULL-d$OFV_RED

#### for-loops ####

loop1 <- lapply(Ninda,function(i) {

loop2 <- lapply(Nindb,function(j) { 

loop3 <- lapply(Nsam,function(n) {

if (i==0) {sama <- 0}
if (i>0) {sama <- sample(d$iOFV[d$FLAG==0],i,replace=TRUE)}                       # Random selection of N ID's from the dataset when FLAG==0

if (j==0) {samb <- 0}
if (j>0) {samb <- sample(d$iOFV[d$FLAG==1],j,replace=TRUE)}                       # Random selection of N ID's from the dataset when FLAG==1

dOFV <- sum(c(sama,samb))

outc <- if(dOFV<=-3.84) {1} else {0}                                            # 1 if alfa is smaller then 0.05 and 0 if bigger then 0.05

})

power <- (length(loop3[loop3==1])/length(loop3))*100

})
})

z <- matrix(unlist(loop1),length(Nindb),length(Ninda))                          # y-power vector to be calculated

### Number of ID's table ###
zzz <- as.data.frame(z)
row.names(zzz) <- Ninda
names(zzz) <- Nindb

wd <- getwd()
write.csv(zzz,paste(wd,"/Cross_results.csv",sep=""))

### Number of ID's 3D plot ###
xx <- rep(Ninda,each=length(Nindb))
yy <- rep(Nindb,length(Ninda))
zz <- as.vector(z)
clc <- rep(NA,length(zz))
clc[zz<power] <- "red"
clc[zz>=power] <- "black"
cl <- rep(1,length(zz))

pdf("plot_crossover.pdf")
s3d <- scatterplot3d(xx,yy,z,
                     type="p",
                     color=rep(clc,cl),
                     highlight.3d=FALSE,
                     angle=55,
                     cex.axis=1.2, 
                     cex.lab=1.2,
                     zlab="power (%)",
                     xlab=designAlab,
                     ylab=designBlab,
                     zlim=c(0,100)
                     )

fit <- c(80,0,0)
s3d$plane3d(fit,col="black",lwd=2)
dev.off()

}