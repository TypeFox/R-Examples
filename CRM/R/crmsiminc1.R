#CRM simulator that called the CRM c function
#codes in this file are used to build CRM program
#This code is used to make CRM package

## function for one simulation; allow one incomplete data 
CRMoneIncomplete <- function(model,nsubject,prior,true,target,a0,b,jump,start.dose,rate,cycle){
  tstart <- rep(NA,nsubject)
  tend <- rep(NA,nsubject)
  lenprior <- as.integer(length(prior))

  # begin trial for the 1st patient
  dose.start <- start.dose   # dose-level at the beginning of trials
  y.rand <- runif(1)
  y <- 1.0*(y.rand<=true[dose.start])  #results for toxicity come from the true distribution
  pData <- c(dose.start,y)
  pData <- matrix(pData,ncol=2)
  i <- 1
  if (y==1) {y.DLT=runif(1,min=0,max=cycle)} # time to DLT
  if (y==0) {y.DLT=cycle}
  t1 <- (rexp(1,rate)) #arrive time 
  tstart[i] <- t1
  tend[i] <- tstart[i] + y.DLT
  new <- cbind(pData,t1,tstart[i],y.DLT,tend[i])
  new <- matrix(new,ncol=6)

  #2nd pt
  modl <- as.integer(model)
  targ <- as.double(target)
  pri <- as.double(prior)
  azero <- as.double(a0)
  bvalue <- as.double(b)
  i <- 2
  callCRM <- .C("CRM",modl,targ,pri,lenprior,azero,bvalue,as.integer(pData),
    as.integer(nrow(pData)),nextDose=as.integer(0),aMean=as.double(0),PACKAGE="CRM")
  x <- callCRM$nextDose
  if (((x-pData[(i-1),1])>1) && (jump==FALSE)) {x<-pData[(i-1),1]+1}
  y.rand <- runif(1)
  y <- 1.0*(y.rand<=true[x])
  pData <- rbind(pData,c(x,y))
  if (y==1) {y.DLT <- runif(1,min=0,max=cycle)}
  if (y==0) {y.DLT <- cycle}
  t1 <- t1+(rexp(1,rate)) #arrive time for pts 2
  if (t1>tend[i-1]) {tstart[i] <- t1}
  if (t1<=tend[i-1]) {tstart[i] <- tend[i-1]}
  tend[i] <- tstart[i] + y.DLT
  new <- rbind(new,c(pData[,1][i],pData[,2][i],t1,tstart[i],y.DLT,tend[i])) 
  
  for (i in 3:nsubject){
    t1 <- t1+(rexp(1,rate))  #arrival time
    #pt minimum start time
    if (t1>tend[i-2]) {tstart[i] <- t1}
    if (t1<=tend[i-2]) {tstart[i] <- tend[i-2]}
    if ((tstart[i])>=tend[i-1]) {
      callCRM <- .C("CRM",modl,targ,pri,lenprior,azero,bvalue,as.integer(pData),
                    as.integer(nrow(pData)),nextDose=as.integer(0),aMean=as.double(0),PACKAGE="CRM")
      x <- callCRM$nextDose
      if (((x-pData[(i-1),1])>1) && (jump==FALSE)) {x <- pData[(i-1),1]+1}
     } else {#if pt came in after end time of last pt then use response of last pt
       tmpData <- matrix(pData[-(i-1),],ncol=2)
       callCRM <- .C("CRM",modl,targ,pri,lenprior,azero,bvalue,as.integer(tmpData),
                     as.integer(nrow(tmpData)),nextDose=as.integer(0),aMean=as.double(0),PACKAGE="CRM")
       x <- callCRM$nextDose
       if (((x-pData[(i-2),1])>1) && (jump==FALSE)) {x <- pData[(i-2),1]+1}
     }
 
    y.rand <- runif(1)
    y <- 1.0*(y.rand<=true[x])
    xy <- cbind(x,y)
    pData <- rbind(pData,xy)
    if (y==1) {y.DLT <- runif(1,min=0,max=cycle)}
    if (y==0) {y.DLT <- cycle}
    tend[i] <- tstart[i] + y.DLT #pt end time
    new=rbind(new,c(pData[,1][i],pData[,2][i],t1,tstart[i],y.DLT,tend[i]))    
  }
  studyTime <- max(tend)
  nTox <- sum(pData[,2])  ## number of toxicity
  doseLevel <- 1:lenprior
  #Proportion of patients treated at each dose level
#  doseTab <- table(pData[,1])/nrow(pData)
  #Proportion of patients treated at each dose level
  doseTab <- table(pData[,1])
  doseTab <- doseTab[match(doseLevel,as.numeric(names(doseTab)))]
  doseTab[is.na(doseTab)] <- 0
  names(doseTab) <- doseLevel
  #Number of toxicities at each dose level
  toxTab <- table(pData[,1],pData[,2])
  if(any(colnames(toxTab) == "1")){
    toxCt <- toxTab[,"1"]
    toxTab <- toxCt[match(doseLevel,as.numeric(rownames(toxTab)))]
    toxTab[is.na(toxTab)] <- 0
  }else {
    toxTab <- rep(0,lenprior)
  }
  names(toxTab) <- doseLevel
  callCRM <- .C("CRM",modl,targ,pri,lenprior,azero,bvalue,as.integer(pData),
                as.integer(nrow(pData)),nextDose=as.integer(0),aMean=as.double(0),PACKAGE="CRM")
  mtd <- callCRM$nextDose
  if (((mtd-pData[nrow(pData),1])>1) && (jump==FALSE)){
    mtd <- pData[nrow(pData),1]+1
  } # make dose escalation <= 1
  result <- list(mtd,nTox,doseTab,toxTab,studyTime)
  return(result)
}

crmsiminc1 <- function(target,prior,true,rate,cycle,nsubject=24,nsim=1000,model=1,a0=1,b=3,jump=FALSE,
                   start.dose=1,seed=777){
  if(target<0 || target > 1){
    stop("Error: target must be greater than 0 and less than 1\n")
  }
  if(any(prior<0) || any(prior>1)){
    stop("Error: All elements in prior must be greater than 0 and less than 1\n")
  }
  if(any(true<0) || any(true>1)){
    stop("Error: All elements in true must be greater than 0 and less than 1\n")
  }     
  if(length(prior)!=length(true)){
    stop("Error: prior and true should have the same length\n")
  }else {
    for(i in 2:length(prior)){
      if(prior[i]<prior[i-1] || true[i]<true[i-1]){
        stop("Error: prior and true must be in an ascending order\n")
      }
    }
  }
  if(model != 1 && model != 2){
    stop("Error: model must be 1 or 2\n")
  }
  if(is.na(match(start.dose,1:length(prior)))){
    stop("Error: start.dose must be in 1:length(prior)\n")
  }
  if(jump != FALSE && jump != TRUE){
    stop("Error: jump must be FALSE or TRUE\n")
  }

  set.seed(seed)
  toxCt <- rep(NA,nsim)
  mtd <- rep(NA,nsim)
  duration <- rep(NA,nsim)
  toxTab <- rep(0,length(prior))
  doseTab <- rep(0,length(prior))
  for (i in 1:nsim){
    fit <- CRMoneIncomplete(model,nsubject,prior,true,target,a0,b,jump,start.dose,rate,cycle)
    mtd[i] <- fit[[1]]
    toxCt[i] <- fit[[2]]
    doseTab <- doseTab + fit[[3]]
    toxTab <- toxTab + fit[[4]]
    duration[i] <- round(fit[[5]])
  }
  doseLevel = 1:length(prior)
  toxTab <- toxTab/nsim
  doseTab <- doseTab/nsim
  propDoseTab <- doseTab/nsubject
  mtdTab <- table(mtd)
  mtdTab <- mtdTab[match(doseLevel,as.numeric(names(mtdTab)))]
  mtdTab[is.na(mtdTab)] <- 0
  mtdTab <- mtdTab/nsim
  names(mtdTab) <- doseLevel
  simTab <- rbind(100*mtdTab,100*propDoseTab,doseTab,toxTab,true)
  simTab <- round(simTab,2)
  rownames(simTab) <- c("% Selection","% Subjects Treated","# Subjects Treated","Average Toxicities",
                        "True probabilities")
  list(SimResult=simTab,TrialDuration=summary(duration))
}

