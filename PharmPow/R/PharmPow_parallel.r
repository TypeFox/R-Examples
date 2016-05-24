PharmPow_parallel <- function(datafullA,
                         datareducedA,
                         datastratifyA,
                         stratifyA,
                         datafullB,
                         datareducedB,
                         datastratifyB,
                         stratifyB,
                         maxIDsA=200,
                         maxIDsB=200,
                         nresamples=10000,
                         power=80,
                         designAlab="Number of patients sampling design A",
                         designBlab="Number of patients sampling design B"){

### read in data

fldns   <- read.table(datafullA,header=TRUE,skip=1)
rddns   <- read.table(datareducedA,header=TRUE,skip=1)
flsps   <- read.table(datafullB,header=TRUE,skip=1)
rdsps   <- read.table(datareducedB,header=TRUE,skip=1)
STRA    <- read.table(datastratifyA,header=TRUE,skip=1)
STRB    <- read.table(datastratifyB,header=TRUE,skip=1)

Ninda   <- seq(0,maxIDsA,by=2)
Nindb   <- seq(0,maxIDsB,by=2)

Nsam    <- seq(1,nresamples,by=1)

### construct sample dataset

dvector1             <- c(fldns$ID,fldns$OBJ,rddns$OBJ)
dA1                  <- matrix(dvector1,length(fldns$ID),3)
dA1                  <- as.data.frame(dA1)
names(dA1)           <- c("ID","OFV_FULL","OFV_RED")
dA1$FLAG             <- 0
dA1$COV              <- 0
STRA2                <- STRA[unique(match(STRA$ID,STRA$ID)),]
dA1$COV              <- STRA2[,paste(stratifyA)]

dvector2             <- c(flsps$ID,flsps$OBJ,rdsps$OBJ)
dA2                  <- matrix(dvector2,length(flsps$ID),3)
dA2                  <- as.data.frame(dA2)
names(dA2)           <- c("ID","OFV_FULL","OFV_RED")
dA2$FLAG             <- 1
dA2$COV              <- 0
STRB2                <- STRB[unique(match(STRB$ID,STRB$ID)),]
dA2$COV              <- STRB2[,paste(stratifyB)]

d                    <- rbind(dA1,dA2)

d$iOFV               <- d$OFV_FULL-d$OFV_RED

#### for-loops

loop1 <- lapply(Ninda,function(i) {

loop2 <- lapply(Nindb,function(j) { 

loop3 <- lapply(Nsam,function(n) {

if (i==0) {sama  <- 0}
if (i>0)  {sama1 <- sample(d$iOFV[d$FLAG==0 & d$COV==0],i/length(unique(d$COV)),replace=TRUE)
           sama2 <- sample(d$iOFV[d$FLAG==0 & d$COV==1],i/length(unique(d$COV)),replace=TRUE)
           sama  <- c(sama1,sama2)}

if (j==0) {samb  <- 0}
if (j>0)  {samb1 <- sample(d$iOFV[d$FLAG==1 & d$COV==0],j/length(unique(d$COV)),replace=TRUE)
           samb2 <- sample(d$iOFV[d$FLAG==1 & d$COV==1],j/length(unique(d$COV)),replace=TRUE)
           samb  <- c(samb1,samb2)}                                             

dOFV <- sum(c(sama,samb))

outc <- if(dOFV<=-3.84) {1} else {0}                                            

})

power <- (length(loop3[loop3==1])/length(loop3))*100

})
})

z <- matrix(unlist(loop1),length(Nindb),length(Ninda))                          

### Number of ID's table
zzz <- as.data.frame(z)
row.names(zzz) <- Ninda
names(zzz) <- Nindb

wd <- getwd()
write.csv(zzz,paste(wd,"/Parallel_results.csv",sep=""))

### Number of ID's 3D plot
xx <- rep(Ninda,each=length(Nindb))
yy <- rep(Nindb,length(Ninda))
zz <- as.vector(z)
clc <- rep(NA,length(zz))
clc[zz<power] <- "red"
clc[zz>=power] <- "black"
cl <- rep(1,length(zz))

pdf("Parallel_results.pdf")
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
