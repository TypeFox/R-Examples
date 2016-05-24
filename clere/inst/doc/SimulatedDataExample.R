source("Methods.R")
## Set seed
Seed  <- 1234
set.seed(Seed)

## Simulation parameters
sigma <-  1
n     <- 25
p     <- 50
g     <-  3

## Beta
intercept <- 0
g         <- 3
probs     <- c(0.36+0.28,0.20,0.12+0.04)
Eff       <- p * probs
a         <- 4
B         <- a**(0:(g-1))-1
Beta      <- rep(B,Eff)

## Generate dataset
nsim  <- 200
N     <- ifelse(n * nsim<5000,5000,n*nsim)
Eps   <- rnorm(N,mean=0,sd=sigma)
xpop  <- matrix(rnorm(N*p),nrow=N,ncol=p)
ypop  <- as.numeric(intercept+xpop%*%Beta+Eps)

numExpSimData <-  NULL
for(isim in 1:nsim){
  cat(paste("\tSimulation #",isim,".\n",sep=""))
  lsim  <- (1+(isim-1)*n):(isim*n)
  xt    <- xpop[+lsim,]; yt <- ypop[+lsim]
  xv    <- xpop[-lsim,]; yv <- ypop[-lsim]
  numExpSimData   <- rbind(numExpSimData,compare(xt,yt,xv,yv,Seed))
}

## Ouptut

save(list="numExpSimData",file="numExpSimData.RData")
save.image("SimulatedDataExample.RData")

## Plot
meths <- c("CLERE0","CLERE","PACS","LASSO",
           "AVG","Ridge","Elastic net",
           "Spike and Slab")
o    <- order(apply(numExpSimData[,1:9],2,median))
dfs  <- round(apply(numExpSimData[,10:18][,o[1:8]],2,mean),1)
sdf  <- round(apply(numExpSimData[,10:18][,o[1:8]],2,sd),1)
tts  <- round(apply(numExpSimData[,19:27][,o[1:8]],2,mean),1)
stt  <- round(apply(numExpSimData[,19:27][,o[1:8]],2,sd),2)
cols <- rainbow(9)

pdf("Simulations.pdf")
par(mar=c(5, 2, 4, 7)+0.1)
boxplot(numExpSimData[,1:9][,o[1:8]],horizontal=TRUE,log="x",
        col=cols,axes=FALSE,pch=18,xlab="Mean Squared Prediction Error")
axis(1)
labs <- paste(meths,"\ndf: ",dfs," (",sdf,")",sep="")
axis(4,at=1:8,labels=labs,las=2)
cts <- paste(tts,"s (",stt,")",sep="")
legend("topleft",legend=cts,box.lty=0,lwd=2,lty=1,col=cols,
       title="Computational time")
dev.off()


