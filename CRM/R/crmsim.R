#CRM simulator that called the CRM c function
#codes in this file are used to build CRM program version 1.2, the latest version
## 8/26/08

##########################################################
CRMsimOne <- function(model,cohort,nsubject,rate,cycle,prior,true,target,a0,b,jump,start.dose){
  lenprior <- as.integer(length(prior))
  dose.start <- rep(start.dose,cohort)   # dose-level at the beginning of trials
  y.rand <- runif(cohort)
  y <- 1.0*(y.rand<=true[dose.start])  #results for toxicity come from the true distribution
  pData <- c(dose.start,y)
  pData <- matrix(pData,ncol=2)
  modl <- as.integer(model)
  targ <- as.double(target)
  pri <- as.double(prior)
  azero <- as.double(a0)
  bvalue <- as.double(b)
  for (i in 2:(nsubject/cohort)){
    callCRM <- .C("CRM",modl,targ,pri,lenprior,azero,bvalue,as.integer(pData),
      as.integer(nrow(pData)),nextDose=as.integer(0),aMean=as.double(0),PACKAGE="CRM")
    x <- callCRM$nextDose
    if (((x-pData[nrow(pData),1])>1) && (jump==FALSE)){
      x <- pData[nrow(pData),1]+1
    } # make dose escalation <= 1
    x.cohort <- rep(x,cohort)
    y.rand <- runif(cohort)
    y <- 1.0*(y.rand<=true[x])
    xy <- cbind(x.cohort,y)
    pData <- rbind(pData,xy)
  }
  nTox<-sum(pData[,2])  ## number of toxicity
  doseLevel <- 1:lenprior
  #Proportion of patients treated at each dose level
#  doseTab <- table(pData[,1])/nrow(pData)
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
  #  mtd <- CRM(model,prior,target,pData,a0,b) #mtd at nsubject+1,proposed mtd for next patient
  callCRM <- .C("CRM",modl,targ,pri,lenprior,azero,bvalue,as.integer(pData),
                as.integer(nrow(pData)),nextDose=as.integer(0),aMean=as.double(0),PACKAGE="CRM")
  mtd <- callCRM$nextDose
  if (((mtd-pData[nrow(pData),1])>1) && (jump==FALSE)){
    mtd <- pData[nrow(pData),1]+1
  } # make dose escalation <= 1

  ### calculate study duration ###
  ## total number of patients 
  ndata <- nrow(pData)
  tarrive <- 0
  tarr <- rep(NA,ndata)
  tstart <- rep(NA,ndata)
  y.DLT <- rep(NA,ndata)
  tend <- rep(NA,ndata)
  
  # for the first cohort of patients
  for (i in 1:cohort){
    if (pData[,2][i]==1) {y.DLT[i] <- runif(1,min=0,max=cycle)}
    if (pData[,2][i]==0) {y.DLT[i] <- cycle}
    #arrival time
    t1 <- (rexp(1,rate)) 
    tarrive <- tarrive+t1
    tarr[i] <- tarrive
    tstart[i] <- tarrive
    tend[i] <- tstart[i] +y.DLT[i]
  }

  #next pts can start after maxend (ie after cohort 1 is done)
  maxend=max(tend[!is.na(tend)])
  #for patients (cohort+1) and on
  j <- 0
  for (i in (cohort+1):ndata){
    j <- j+1
    t1 <- (rexp(1,rate)) 
    tarrive <- tarrive+t1 #arrival time
    tarr[i] <- tarrive
    if (tarr[i]>maxend) {tstart[i] <- tarrive}  #find start time
    if (tarr[i]<=maxend) {tstart[i] <- maxend}
    if (pData[,2][i]==1) {y.DLT[i] <- runif(1,min=0,max=cycle)} #time to DLT
    if (pData[,2][i]==0) {y.DLT[i] <- cycle}
    tend[i] <- tstart[i] +y.DLT[i]  #end time of each pt
    if (j==cohort){
      maxend <- max(tend[!is.na(tend)])
      j <- 0
    }
  }
#  new=cbind(tarr,tstart,y.DLT,tend)
#  print(new)
  studydur <- max(tend)
  result <- list(mtd,nTox,doseTab,toxTab,studydur)
  return(result)
}
#####################################

crmsim <- function(target,prior,true,rate,cycle,cohort=1,nsubject=24,nsim=1000,model=1,a0=1,b=3,jump=FALSE,
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
  if(rate<0 || cycle <0){
    stop("Error: rate and cycle must be greater than 0\n")
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
  if(floor(nsubject/cohort) != (nsubject/cohort)){
    mess = paste("nsubject/cohort is not an integer. The number of subjects is changed to",
      cohort*floor(nsubject/cohort))
    warning(mess);
  }
  set.seed(seed)
  toxCt <- rep(NA,nsim)
  mtd <- rep(NA,nsim)
  duration <- rep(NA,nsim)
  toxTab <- rep(0,length(prior))
  doseTab <- rep(0,length(prior))
  for (i in 1:nsim){
    fit <- CRMsimOne(model,cohort,nsubject,rate,cycle,prior,true,target,a0,b,jump,start.dose)
    mtd[i] <- fit[[1]]
    toxCt[i] <- fit[[2]]
    doseTab <- doseTab + fit[[3]]
    toxTab <- toxTab + fit[[4]]
    duration[i] <- round(fit[[5]])
  }
  doseLevel <- 1:length(prior)
  toxTab <- toxTab/nsim
  doseTab <- doseTab/nsim
  propDoseTab <- doseTab/(cohort*floor(nsubject/cohort))
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

