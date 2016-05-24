library("MASS")
library("parallel")
library("clere") 

Seed      <-  1234 
set.seed(Seed)

nsim      <-   200 # set nsim > 1
plotit    <- FALSE
sparse    <- FALSE
nItEM     <-  2000
nsamp     <-  5000
analysis  <-  "fit"
maxit     <- 10000
tol       <-  1e-8
n         <-    25
p         <-    50
intercept <-     0
sigma     <-     1
gamma     <-     0

g     <- 3
probs <- c(0.36+0.28,0.20,0.12+0.04)
Eff   <- p * probs
a     <- 4
B     <- a**(0:(g-1))-1
Beta  <- rep(B,Eff)
u     <- rep(NA,nsim)

Pred  <- cbind(MCEMA=u,MCEMB=u,MCEMC=u,SEM=u,ORACLE=u)
Time  <- cbind(MCEMA=u,MCEMB=u,MCEMC=u,SEM=u)
Bias  <- cbind(MCEMA=u,MCEMB=u,MCEMC=u,SEM=u)
Liks  <- cbind(MCEMA=u,MCEMB=u,MCEMC=u,SEM=u,ORACLE=u)


mA <-   5
mB <-  25
mC <- 125

perms <- cbind(c(1,2,3),
               c(1,3,2),
               c(2,1,3),
               c(2,3,1),
               c(3,1,2),
               c(3,2,1))

thetaTrue <- c(intercept,B,probs,sigma^2,gamma^2)
ThetaA    <- matrix(NA,nrow=nsim,ncol=length(thetaTrue))
colnames(ThetaA) <- c("beta0",paste("b",1:g,sep=""),paste("pi",1:g,sep=""),"sigma2","gamma2")
ThetaB    <- ThetaA
ThetaC    <- ThetaA
ThetaS    <- ThetaA


gtheta <- function(theta){
  inds <- c(1,8,9)
  subtheta <- theta[-inds]
  j <- which.min(apply(perms,2,function(j) sum((thetaTrue[-inds]-c(subtheta[1:3][j],subtheta[4:6][j]))^2) ))
  return(c(theta[1],theta[2:4][perms[,j]],theta[5:7][perms[,j]],theta[8:9]))
}

## Prediction
N     <- ifelse(n * nsim<1000,1000,n*nsim)
Eps   <- rnorm(N,mean=0,sd=sigma)
xpop  <- matrix(rnorm(N*p),nrow=N,ncol=p)
ypop  <- as.numeric(intercept+xpop%*%Beta+Eps)

## Parallel 
nCPU     <- 10 ## was run with nCPU = 10 in the submitted article
parallel <- TRUE

for(sim in 1:nsim){
  
  cat(paste("\tSimulation n° ",sim,".\n",sep=""))
  lsim  <- (1+(sim-1)*n):(sim*n)
  X     <- xpop[+lsim,]; Y  <- ypop[+lsim]
  Xv    <- xpop[-lsim,]; Yv <- ypop[-lsim]
  
  ## MCEM A
  mcttA  <- system.time ( mcmodA <- fitClere(y=Y,x=X,g=g,
              analysis=analysis,parallel=parallel,nstart=nCPU,
              algorithm="MCEM",seed=Seed,
              plotit=plotit,sparse=sparse,
              nItEM=nItEM,nBurn=max(1,0.2*mA),
              nItMC=1,nsamp=mA) )
  mcthetaA  <- gtheta( c(mcmodA@intercept,mcmodA@b,mcmodA@pi,mcmodA@sigma2,mcmodA@gamma2) )
  zstartA   <- clusters(mcmodA)
  misgrpA   <- which(!(1:g)%in%zstartA)
  if(length(misgrpA)){
    tabA <- table(zstartA)
    zstartA[which(zstartA==as.numeric(names(tabA[which.max(tabA)])))[1:length(misgrpA)]] <- misgrpA
  }
  tmpModA   <- fitClere(y=Y,x=X,g=g,analysis=analysis,parallel=parallel,nstart=nCPU,
                        algorithm="SEM",plotit=plotit,sparse=sparse,
                        nItEM=1,nBurn=0,nItMC=0,nsamp=nsamp,seed=Seed,
                        theta0=mcthetaA,Z0=zstartA-1)
  ## MCEM B
  mcttB  <- system.time ( mcmodB <- fitClere(y=Y,x=X,g=g,
              analysis=analysis,parallel=parallel,nstart=nCPU,
              algorithm="MCEM",seed=Seed,
              plotit=plotit,sparse=sparse,
              nItEM=nItEM,nBurn=max(1,0.2*mB),
              nItMC=1,
              nsamp=mB))  
  mcthetaB  <- gtheta( c(mcmodB@intercept,mcmodB@b,mcmodB@pi,mcmodB@sigma2,mcmodB@gamma2) )
  zstartB   <- clusters(mcmodB)
  misgrpB   <- which(!(1:g)%in%zstartB)
  if(length(misgrpB)){
    tabB <- table(zstartB)
    zstartB[which(zstartB==as.numeric(names(tabB[which.max(tabB)])))[1:length(misgrpB)]] <- misgrpB
  }
  tmpModB   <- fitClere(y=Y,x=X,g=g,analysis=analysis,parallel=parallel,nstart=nCPU,
                        algorithm="SEM",plotit=plotit,sparse=sparse,
                        nItEM=1,nBurn=0,nItMC=0,nsamp=nsamp,seed=Seed,
                        theta0=mcthetaB,Z0=zstartB-1)
  
  ## MCEM C
  mcttC  <- system.time ( mcmodC <- fitClere(y=Y,x=X,g=g,
              analysis=analysis,parallel=parallel,nstart=nCPU,
              algorithm="MCEM",seed=Seed,
              plotit=plotit,sparse=sparse,
              nItEM=nItEM,nBurn=max(1,0.2*mC),
              nItMC=1,
              nsamp=mC))
  mcthetaC  <- gtheta( c(mcmodC@intercept,mcmodC@b,mcmodC@pi,mcmodC@sigma2,mcmodC@gamma2) )
  zstartC   <- clusters(mcmodC)
  misgrpC   <- which(!(1:g)%in%zstartC)
  if(length(misgrpC)){
    tabC <- table(zstartC)
    zstartC[which(zstartC==as.numeric(names(tabC[which.max(tabC)])))[1:length(misgrpC)]] <- misgrpC
  }
  tmpModC   <- fitClere(y=Y,x=X,g=g,analysis=analysis,parallel=parallel,nstart=nCPU,
                        algorithm="SEM",plotit=plotit,sparse=sparse,
                        nItEM=1,nBurn=0,nItMC=0,nsamp=nsamp,seed=Seed,
                        theta0=mcthetaC,Z0=zstartC-1)
    
  ## SEM 
  tt2  <- system.time ( mod2 <- fitClere(y=Y,x=X,g=g,
                        analysis=analysis,parallel=parallel,nstart=nCPU,
                        algorithm="SEM",seed=Seed,
                        plotit=plotit,sparse=sparse,
                        nItEM=nItEM,nBurn=nItEM/2,
                        nItMC=10,nsamp=nsamp))
  theta2 <- gtheta( c(mod2@intercept,mod2@b,mod2@pi,mod2@sigma2,mod2@gamma2) )

  ## oracle
  thetaOracle <- thetaTrue
  thetaOracle[9] <- 1e-16
  zOracle <- rep(0:2,Eff)
  oracle  <- fitClere(y=Y,x=X,g=g,analysis=analysis,parallel=parallel,nstart=nCPU,
                       algorithm="SEM",plotit=plotit,sparse=sparse,
                       nItEM=1,nBurn=0,nItMC=0,nsamp=nsamp,seed=Seed,
                       theta0=thetaOracle,Z0=zOracle)
  
  Pred[sim,"MCEMA"]   <- mean((Yv-predict(mcmodA ,Xv))^2,na.rm=TRUE)
  Pred[sim,"MCEMB"]   <- mean((Yv-predict(mcmodB ,Xv))^2,na.rm=TRUE)
  Pred[sim,"MCEMC"]   <- mean((Yv-predict(mcmodC ,Xv))^2,na.rm=TRUE)
  Pred[sim,"SEM"]     <- mean((Yv-predict(mod2   ,Xv))^2,na.rm=TRUE)
  Pred[sim,"ORACLE"]  <- mean((Yv-predict(oracle ,Xv))^2,na.rm=TRUE)
  
  Time[sim,"MCEMA"]   <- mcttA["elapsed"]
  Time[sim,"MCEMB"]   <- mcttB["elapsed"]
  Time[sim,"MCEMC"]   <- mcttC["elapsed"]
  Time[sim,"SEM"]     <- tt2["elapsed"]
  
  Bias[sim,"MCEMA"]  <- sum((thetaTrue-mcthetaA)^2)
  Bias[sim,"MCEMB"]  <- sum((thetaTrue-mcthetaB)^2)
  Bias[sim,"MCEMC"]  <- sum((thetaTrue-mcthetaC)^2)
  Bias[sim,"SEM"]    <- sum((thetaTrue-theta2)^2)
  
  Liks[sim,"MCEMA"]  <- tmpModA@likelihood
  Liks[sim,"MCEMB"]  <- tmpModB@likelihood
  Liks[sim,"MCEMC"]  <- tmpModC@likelihood
  Liks[sim,"SEM"]    <- mod2@likelihood
  Liks[sim,"ORACLE"] <- oracle@likelihood

  ThetaA[sim,] <- mcthetaA
  ThetaB[sim,] <- mcthetaB
  ThetaC[sim,] <- mcthetaC
  ThetaS[sim,] <- theta2
}

## Output
algoComp <- list(Bias=Bias,Pred=Pred,Time=Time,Liks=Liks)
save(list="algoComp",file="algoComp.RData")
save.image("SEM_vs_MCEM_simulations.RData")

## Latex table
createTable1 <- function(algoComp){
  nsim   <- 200
  lnames <- names(algoComp)
  for(i in 1:length(lnames)){
    Tmp <- algoComp[[i]][,1:4]
    colnames(Tmp) <- c("MCEM5","MCEM25","MCEM125","SEM")
    if(lnames[i]=="Liks"){
      tmp <- format(100*table(factor(apply(Tmp,1,which.max),levels=1:4))/nrow(Tmp),digit=2,nsmall=1)
      names(tmp) <- colnames(Tmp)
    }else{
      tmp <- 100*table(factor(apply(Tmp,1,which.min),levels=1:4))/nrow(Tmp)
      names(tmp) <- colnames(Tmp)
    }
    assign(paste("isBest",lnames[i],sep=""),tmp)
    ## Mean
    avTmp <- format(apply(Tmp,2,median),digit=2,nsmall=1)
    sdTmp <- format(1.253*apply(Tmp,2,sd)/sqrt(nsim),digit=2,nsmall=1)
    assign(paste("av",lnames[i],sep=""),avTmp)
    assign(paste("sd",lnames[i],sep=""),sdTmp)
  }
  avTPpred <- format( median(algoComp$Pred[,"ORACLE"]),digit=2,nsmall=1)
  sdTPpred <- format( 1.253*sd(algoComp$Pred[,"ORACLE"])/sqrt(nsim),digit=2,nsmall=1)
  avTPliks <- format( median(algoComp$Liks[,"ORACLE"]),digit=2,nsmall=1)
  sdTPliks <- format( 1.253*sd(algoComp$Liks[,"ORACLE"])/sqrt(nsim),digit=2,nsmall=1)
  
  tab <- paste(
               "\\begin{center}\n",
               "\\begin{table}[h!]\n",
               "\\begin{tabular}{llrr}",
               "\\toprule\n",
               "              &                  & \\% of times                    & Median       \\\\\n",
               "Performance indicators & Algorithms & the algorithm was best       & (Std. Err.)\\\\\n",
               "\\midrule\n",
               "CT (seconds)  &  SEM             &  ",isBestTime["SEM"],"           &  ",avTime["SEM"]," ( ",sdTime["SEM"]," ) \\\\\n",
               "              &  MCEM$_5$        &  ",isBestTime["MCEM5"],"         &  ",avTime["MCEM5"]," ( ",sdTime["MCEM5"]," ) \\\\\n",
               "              &  MCEM$_{25}$     &  ",isBestTime["MCEM25"],"        &  ",avTime["MCEM25"]," ( ",sdTime["MCEM25"]," ) \\\\\n",
               "              &  MCEM$_{125}$    &  ",isBestTime["MCEM125"],"       &  ",avTime["MCEM125"]," ( ",sdTime["MCEM125"]," ) \\\\\n",
               "\\\\\n",
               "\\midrule\n",
               "MSEE          &  SEM             &  ",isBestBias["SEM"],"           &  ",avBias["SEM"]," ( ",sdBias["SEM"]," ) \\\\\n",
               "              &  MCEM$_5$        &  ",isBestBias["MCEM5"],"         &  ",avBias["MCEM5"]," ( ",sdBias["MCEM5"]," ) \\\\\n",
               "              &  MCEM$_{25}$     &  ",isBestBias["MCEM25"],"        &  ",avBias["MCEM25"]," ( ",sdBias["MCEM25"]," ) \\\\\n",
               "              &  MCEM$_{125}$    &  ",isBestBias["MCEM125"],"       &  ",avBias["MCEM125"]," ( ",sdBias["MCEM125"]," ) \\\\\n",
               "\\\\\n",
               "\\midrule\n",
               "MSPE          &  SEM             &  ",isBestPred["SEM"],"           &  ",avPred["SEM"]," ( ",sdPred["SEM"]," ) \\\\\n",
               "              &  MCEM$_5$        &  ",isBestPred["MCEM5"],"         &  ",avPred["MCEM5"]," ( ",sdPred["MCEM5"]," ) \\\\\n",
               "              &  MCEM$_{25}$     &  ",isBestPred["MCEM25"],"        &  ",avPred["MCEM25"]," ( ",sdPred["MCEM25"]," ) \\\\\n",
               "              &  MCEM$_{125}$    &  ",isBestPred["MCEM125"],"       &  ",avPred["MCEM125"]," ( ",sdPred["MCEM125"]," ) \\\\\n",                             
               "              &  True parameter  &  ---                             &  ",avTPpred," (",sdTPpred," )\\\\\n",
               "\\\\\n",
               "\\midrule\n",
               "ML            &  SEM             &  ",isBestLiks["SEM"],"           &  ",avLiks["SEM"]," ( ",sdLiks["SEM"]," ) \\\\\n",
               "              &  MCEM$_5$        &  ",isBestLiks["MCEM5"],"         &  ",avLiks["MCEM5"]," ( ",sdLiks["MCEM5"]," ) \\\\\n",
               "              &  MCEM$_{25}$     &  ",isBestLiks["MCEM25"],"        &  ",avLiks["MCEM25"]," ( ",sdLiks["MCEM25"]," ) \\\\\n",
               "              &  MCEM$_{125}$    &  ",isBestLiks["MCEM125"],"       &  ",avLiks["MCEM125"]," ( ",sdLiks["MCEM125"]," ) \\\\\n",        
               "              &  True parameter  &  ---                             &  ",avTPliks," (",sdTPliks," )\\\\\n",
               "\\bottomrule\n",
               "\\end{tabular}\n",
               "\\caption{\\label{tab:simulations} Performance indicators used to compare SEM and MCEM algorithms. Computational Time (CT) was measured on a Intel(R) Xeon(R) CPU E7- 4870  @ 2.40GHz processor. The best algorithm is defined as the one that either reached the largest log-likelihood (ML) or the lowest CT, Mean Squared Prediction Error (MSPE) and Mean Squared Estimation Error (MSEE).}\n",
               "\\end{table}\n",
               "\\end{center}\n")
  return(tab)  
}

cat( createTable1(algoComp) )

