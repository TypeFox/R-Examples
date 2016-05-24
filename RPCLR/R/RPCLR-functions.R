#### Main function is called "GetVarImp()
#### Inputs to this function are:
#### (1) MyData: Dataset with n (number of subjects) rows and p (number of features) columns
#### (2) MyOut: Vector of length n (number of subjects) containing binary indicators of outcome (Case vs control) for each subject
#### (3) MyStrat: Vector of length n (number of subjects) containing matched stratum indicator
#### (4) mtry: Number of variables (Range: 1 to p) to be randomly sampled for inclusion in each penalized model
#### (5) numBS: Number of bootstrap datasets to be sampled
#### Output of this function: Vector of importance scores for each feature.

#### The main function GetVarImp() calls the following sub-functions:
#### (1) CreateBSDataset() - obtains a bootstrap dataset of the same size as the original dataset, by sampling the matched pairs with replacement
#### (2) GetIntData() - generates all pairwise interactions between the set of features input to the function
#### (3) FitPenalizedModel() - fits the penalized conditional logistic regression model
#### (4) GetOOBAIC() - gets the AIC based on data from the 'out of bag (OOB)' subjects

#### The function GenerateData() generates a simulated dataset of size n (number of subjects ) x p (number of features).
#### Inputs to this function are:
#### (1) numstrat - number of matched pairs
#### (2) NumType.BM - number of biomarkers (with non-zero difference in means between cases and controls)
#### (3) NumType.NS - number of noise features (with no mean difference between cases and controls)
#### (4) mu.diff - difference in means between cases and controls for biomarkers
#### (5) rho - correlation between matched pair (for biomarkers only)
#### Output of this function is a list with elements: 
#### (1) Data: n x p data matrix
#### (2) Out: n x 1 vector of outcomes (Case versus control)
#### (3) Strat: n x 1 vector of matched stratum indicators 

GetVarImp <- function(MyData, MyOut, MyStrat, mtry,numBS){

CreateBSDataset <- function(MyData, MyOut, MyStrat, index){
  NewData <- c()
  NewOut <- c()
  NewStrat <- c()
  for (i in 1:length(index)){
    ind <- which(MyStrat == index[i])
    NewData <- rbind(NewData, MyData[ind,])
    NewOut <- c(NewOut, MyOut[ind])
    NewStrat <- c(NewStrat,rep(i,2)) 
  }

  out <- list(NewData=NewData, NewOut=NewOut,NewStrat=NewStrat)
  return(out)
}


GetIntData<- function(BSData){
  IntData <- c()
  mynames <- colnames(BSData)
  intnames <- c()
   for (i in 2:ncol(BSData)){
    for (j in 1:(i-1)){
      IntData <- cbind(IntData, BSData[,i]*BSData[,j])
      intnames <- as.character(c(intnames, paste(mynames[i],"-",mynames[j],sep="")))
          }
  }
  colnames(IntData) <- intnames
  return(IntData)
}


FitPenalizedModel <- function(Data,Out, Strat){
        numvar <- ncol(Data)
    	MyMod <- clogit(Out ~ ridge(Data[,1:numvar])+ strata(Strat)) 		
     	mybeta <-   MyMod$coef    	
  	return(mybeta)
}

GetOOBAIC <- function(mymod.beta, MyData, MyOut, MyStrat, index, covind, createInt=TRUE){
  allstrat <- sort(unique(MyStrat))
  instrat <- sort(unique(index))
  outstrat <- allstrat[-instrat]
  coefs <- as.numeric(mymod.beta)

  ### Create data of pairwise interactions for selected covariates
  Fun.ModData <- MyData[,covind]
  if (createInt){
    Fun.IntData <- GetIntData(MyData[,covind])
    Fun.ModData <- cbind(MyData[,covind], Fun.IntData)}
  OOB.Lik <- 0
 
  for (i in 1:length(outstrat)){
    this1 <- Fun.ModData[MyStrat == outstrat[i] & MyOut == 1,]
    this0 <- Fun.ModData[MyStrat == outstrat[i] & MyOut == 0,]
    #print(length(coefs))
    #print(length(this1))
    pred1 <- exp(coefs%*%this1)
    pred0 <- exp(coefs%*%this0)
    outi <-  pred1/(pred0 + pred1)   
    OOB.Lik <- OOB.Lik + log(outi) }

  out <- OOB.Lik
  return(out)
}


  ##### main body of function
  library(survival)
  
  numstrat <- length(unique(MyStrat))
  out <- matrix(NA, nrow=numBS, ncol=ncol(MyData))  
  for (i in 1:numBS){
     ### Obtain bootstrap dataset
      index <- sample(1:numstrat, numstrat, replace=TRUE)
      Out <- CreateBSDataset(MyData, MyOut, MyStrat, index) 
      NewData <- Out$NewData
      NewOut <- Out$NewOut
      NewStrat <- Out$NewStrat
      fts.names <- colnames(NewData)
      ### Obtain random subset of covariates 
      covind <- sample(1:ncol(NewData), mtry, replace=FALSE)
      BSData <- NewData[,covind]
      BSOut <- NewOut
      BSStrat <- NewStrat

     ### Augment BS data to get all pairwise interactions
    IntData <- GetIntData(BSData)
    ModData <- cbind(BSData, IntData)

    #####################################################################     
    ### Fit penalized CLR models using ridge penalty
    #### Use default option for ridge penalty
    ######################################################################
	print(paste("BS-", i, sep=""))

    mymod.beta <- FitPenalizedModel(ModData, BSOut, BSStrat)   
    fullAIC <- GetOOBAIC(mymod.beta, MyData, MyOut, MyStrat, index, covind)
    subAIC <- c()
    for (j in 1:mtry){   			
		IntDataj <- GetIntData(BSData[,-j])  	
		ModDataj <- cbind(BSData[,-j], IntDataj)  
		mymodS <- FitPenalizedModel(ModDataj, BSOut, BSStrat)
		AICj <- GetOOBAIC(mymodS, MyData, MyOut, MyStrat, index, covind[-j])   				    
		subAIC <- c(subAIC, AICj)}
		
    ChangeLik <- (fullAIC-subAIC)
    out[i,covind] <- ChangeLik
  }
  out[out == "NaN"] <- NA
  myfun <- function(vec){
  	vec[vec == Inf] <- NA
  	vec[vec == -Inf] <- NA
  	return(mean(vec[!is.na(vec)]))}
  meanout <- apply(out, 2, myfun) 
  return(meanout)  
}


GenerateData <- function(numstrat, NumType.BM, NumType.NS, mu.diff, rho){
  library(MASS)
  MyData <- c()
  Datai.BM <- c()
  Datai.NS <- c()

  cov.mat <- rbind(c(1,rho), c(rho,1))
  for (j in 1:NumType.BM){
    temp <- mvrnorm(n=numstrat, mu=c(1,1+mu.diff), Sigma=cov.mat)
    newtemp <- c(temp[,1], temp[,2])
    Datai.BM <- cbind(Datai.BM, newtemp)
  } 
  My.Colnames1 <- rep(paste("E-", 1+mu.diff,"-C-",rho, sep=""), NumType.BM)

  cov.mat2 <- rbind(c(1,0), c(0,1))
  for (j in 1:NumType.NS){
    temp <- mvrnorm(n=numstrat, mu=c(1,1), Sigma=cov.mat2)
    newtemp <- c(temp[,1], temp[,2])
    Datai.NS <- cbind(Datai.NS, newtemp)
  } 
  My.Colnames2 <- rep(paste("E-", 1,"-C-",0, sep=""), NumType.NS)

  MyData <- cbind(Datai.BM, Datai.NS) 
  colnames(MyData) <- c(as.character(My.Colnames1), as.character(My.Colnames2))
  MyOut <- c(rep(0,numstrat), rep(1,numstrat))
  MyStrat <- c(seq(1,numstrat), seq(1,numstrat))

  out <- list(Data=MyData, Out = MyOut, Strat=MyStrat)
  return(out) 
}


