## Library and scripts
library(parallel)
library(glmnet)
library(lars)
library(spikeslab)
library(clere)

## Implementation of the AVG methodology
collapse <- function(x,classes){
  h <- t(apply(x,1,function(u) aggregate(x=u,by=list(classes),sum)[,"x"] ))
  return(h)
}

avg <- function(x,y,nCPU){
  ## hierarchical clustering
  p   <- ncol(x)
  don <- scale(x, center = TRUE, scale = TRUE)
  dc  <- dist(t(don), method ="euclidean", diag=FALSE, upper=FALSE)
  hier<- hclust(dc,"ward.D")
  sgs <- lapply(1:p, function(j) cutree(hier,j))[-1]
  bhi <- mclapply(sgs,function(classes) min( cv.glmnet(collapse(x,classes),as.vector(y),type="mse",alpha=1.0,nfolds=5)$cvm ),
                  mc.cores=nCPU,mc.set.seed=TRUE)
  classes <- sgs[[which.min(bhi)]]
  g <- length(unique(classes))
  z <- matrix(0,nrow=p,ncol=g)
  for(j in 1:p){
    z[j,classes[j]] <- 1
  }
  xz <- x%*%z
  cv <-  cv.glmnet(xz,as.vector(y),type="mse",alpha=1.0,nfolds=5)
  lminlasso <- cv$lambda.min
  fitlasso  <- glmnet(xz,as.vector(y),lambda=lminlasso,alpha=1.0,standardize=FALSE)
  Bavg      <- z%*%fitlasso$beta[,1]
  df        <- length(unique(Bavg[which(Bavg!=0)]))
  return( list(Bavg = Bavg, Z = z, a0 = fitlasso$a0, df = df) )
}

compare <- function(xt,yt,xv,yv,seed=seed,nCPU=10){
  set.seed(seed)

  Competitors  <- c("LASSO","RIDGE","ELNET","STEP","CLERE","CLERE_s","SS","AVG","PACS",
                    "plasso","pridge","pelnet","pstep","pclere","pclere_s","pss","pavg","ppacs",
                    "tlasso","tridge","telnet","tstep","tclere","tclere_s","tss","tavg","tpacs",
                    "seed")
  nComp          <- length(Competitors)
  Pred           <- matrix(NA,nrow=1,ncol=nComp)
  colnames(Pred) <- Competitors
  sim <- 1
  Pred[sim,"seed"] <- seed
  ## STEP
  tstep <- system.time({
                       stepwise <- lars(xt,yt,type="stepwise",intercept=TRUE,normalize=FALSE);                     
                       B        <- stepwise$beta}
                       )
  Bstep = B[nrow(B),]
  Pred[sim,"STEP"]  <- mean((yv-stepwise$mu-xv%*%Bstep)^2,na.rm=TRUE)
  Pred[sim,"pstep"] <- sum(Bstep!=0)
  Pred[sim,"tstep"] <- tstep[3]
  
  ## LASSO
  tlasso <- system.time({
    cvlasso   <- cv.glmnet(xt,as.vector(yt),type="mse",alpha=1.0,nfolds=5)
    lminlasso <- cvlasso$lambda.min
    fitlasso  <- glmnet(xt,as.vector(yt),lambda=lminlasso,alpha=1.0)
    Blasso    <- fitlasso$beta[,1]
  })
  Pred[sim,"LASSO"]  <- mean((yv-fitlasso$a0-xv%*%Blasso)^2,na.rm=TRUE)
  Pred[sim,"plasso"] <- sum(Blasso!=0)
  Pred[sim,"tlasso"] <- tlasso[3]
  
  ## RIDGE    
  tridge <- system.time({
    cv        <- cv.glmnet(xt,as.vector(yt),type="mse",alpha=0.0,nfolds=5)
    lminridge <- cv$lambda.min
    fitridge  <- glmnet(xt,as.vector(yt),lambda=lminridge,alpha=0.0)
    Bridge    <- fitridge$beta[,1]
  })
  Pred[sim,"RIDGE"]   <- mean((yv-fitridge$a0-xv%*%Bridge)^2,na.rm=TRUE)
  Pred[sim,"pridge"]  <- sum(Bridge!=0)
  Pred[sim,"tridge"]  <- tridge[3]
  
  ## ELNET
  telnet <- system.time({
    aset  <- seq(0.0,0.0,by=0.1)
    laset <- length(aset)
    cven  <- matrix(NA,nrow=laset,ncol=2)
    for(k in 1:laset){
      cv <- cv.glmnet(xt,as.vector(yt),type="mse",alpha=aset[k],nfolds=5)
      cven[k,1] <- min(cv$cvm)
      cven[k,2] <- cv$lambda.min
    }
    imin <- which.min(cven[,1])
    fitelnet  <- glmnet(xt,as.vector(yt),lambda=cven[imin,2],alpha=aset[imin])
    Belnet    <- fitelnet$beta[,1]
  })
  Pred[sim,"ELNET"]  <- mean((yv-fitelnet$a0-xv%*%Belnet)^2,na.rm=TRUE)
  Pred[sim,"pelnet"] <- sum(Belnet!=0)
  Pred[sim,"telnet"] <- telnet[3]

  ## AVG
  ttavg            <- system.time( avgmod <- avg(xt,yt,nCPU) )
  Pred[sim,"AVG"]  <- mean((yv-avgmod$a0-xv%*%avgmod$Bavg)^2,na.rm=TRUE)
  Pred[sim,"pavg"] <- avgmod$df
  Pred[sim,"tavg"] <- ttavg[3]  
  
  ## CLERE
  gmax   <- 3
  nstart <- nCPU
  w    <- function(gmax,sparse){
    mod <- fitClere(y=yt,x=xt,g=gmax,seed=seed,
                    nstart=nstart,parallel=TRUE,
                    analysis="aic",plotit=FALSE,
                    sparse=sparse,nItEM=2000,
                    nBurn=1000,nItMC=10,dp=5,
                    nsamp=1000)
    return(mod)
  }
  tclere0 <- system.time( clere0 <- w(gmax,FALSE) )
  tclere1 <- system.time( clere1 <- w(gmax,TRUE)  )
  
  Pred[sim,"CLERE"]    <- mean((yv-predict(clere0,xv))^2,na.rm=TRUE)
  Pred[sim,"pclere"]   <- 2*length(clere0@b)+2
  Pred[sim,"tclere"]   <- tclere0[3]
  
  Pred[sim,"CLERE_s"]  <- mean((yv-predict(clere1,xv))^2,na.rm=TRUE)
  Pred[sim,"pclere_s"] <- 2*length(clere1@b)+1
  Pred[sim,"tclere_s"] <- tclere1[3]

  ## Spike and Slab
  tss <- system.time( ss  <- spikeslab(x = xt, y = yt,n.iter1 = 1000, n.iter2 = 1000, seed = as.integer(seed)) )
  pss <- predict.spikeslab(ss,newdata=xv)    
  Pred[sim,"SS"]  <- mean((yv-pss$yhat.bma)^2)
  Pred[sim,"pss"] <- ss$phat
  Pred[sim,"tss"] <- tss[3]

  ## PACS
  epsPACS <- 1e-5  
  lpacs   <- c(outer(c(1,2,5),10**(-2:2),FUN="*"))
  nfold   <- 5
  dp      <- round(nrow(xt)/nfold)

  pp <- function(lambda){
    fitpacs     <- fitPacs(yt,xt,lambda,betaInput=Bridge,epsPACS=epsPACS,nItMax=1000)
    return(fitpacs)
  }
  ih <- function(ifold,lambda){
    lfold        <- (1+(ifold-1)*dp):(min(ifold*dp,nrow(xt)))
    fitpacs      <- fitPacs(yt[-lfold],xt[-lfold,],lambda,betaInput=Bridge,epsPACS=epsPACS,nItMax=1000)
    return(mean((yt[lfold]-fitpacs@a0-xt[lfold,]%*%fitpacs@betaOutput)^2))
  }
  
  cvpp <- function(lambda){
    cvout <- rep(NA,nfold)
    for(ifold in 1:nfold){
      cvout[ifold] <- ih(ifold,lambda)
    }
    return(sqrt(mean(as.numeric(cvout),na.rm=TRUE)))
  }
  ttpacs <- system.time({
                        cvpacs    <- do.call("c",mclapply(lpacs,cvpp,mc.cores=length(lpacs),mc.set.seed=TRUE))
                        lpacsGood <- lpacs[which.min(cvpacs)[1]]
                        pacsGood  <- pp(lpacsGood)
                        })
  
  Pred[sim,"PACS"]  <- mean((yv-pacsGood@a0-xv%*%pacsGood@betaOutput)^2,na.rm=TRUE)
  Pred[sim,"ppacs"] <- pacsGood@K
  Pred[sim,"tpacs"] <- ttpacs[3]
  
  return(Pred)
}
