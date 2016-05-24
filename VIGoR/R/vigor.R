vigor<-function(Pheno,Geno,Method=c("BL","EBL","wBSR","BayesB","BayesC","SSVS","MIX"),Hyperparameters,
                Function="fitting",Nfold=10,CVFoldTuning=5,Partition=NULL,
                Covariates="Intercept",Threshold=2+log10(ncol(Geno)),Maxiterations=1000,RandomIni=FALSE,Printinfo=TRUE){

  #functions
  genomewideregression<-function(Pheno,Geno,Methodcode,Hyperparameters,Covariates,Threshold,Maxiterations,Randomize){
    
    N<-length(Pheno)
    P<-ncol(Geno) 
    if(Methodcode==1|Methodcode==2){Priortype<-1} else {Priortype<-2}
    Fc<-ncol(Covariates)  
    Threshold<-10^(-Threshold)
    
    CondResidual<-1
    expDelta2<-0
    Tau0<-numeric(2)
    LBmonitor<-Rmonitor<-numeric(Maxiterations)
    Expectations<-Uncertainty<-Variance<-Gamma<-Eta2<-numeric(P)
    ExpectationsQ<-UncertaintyQ<-numeric(Fc)
    
    Result<-.C("vigorr",
               as.integer(Priortype),
               as.integer(Methodcode),
               as.integer(CondResidual),
               as.integer(P),
               as.integer(Fc),
               as.integer(N),
               as.integer(Maxiterations),
               as.integer(Randomize),
               as.double(Threshold),
               as.double(Pheno),
               as.double(Geno[1:(length(Geno))]),
               as.double(Covariates[1:length(Covariates)]),
               as.double(Hyperparameters),
               as.double(Tau0),
               as.double(LBmonitor),
               as.double(Rmonitor),
               as.double(Expectations),
               as.double(Uncertainty),
               as.double(Variance),
               as.double(ExpectationsQ),
               as.double(UncertaintyQ),
               as.double(Gamma),
               as.double(expDelta2),            
               as.double(Eta2),
               PACKAGE="VIGoR"
    ) 
    
    for(i in 1:14){Result[1]<-NULL};
    if(Methodcode==1){
      Result[8]<-NULL; Result[9]<-NULL;
      names(Result)<-c("LB","ResidualVar","Beta","Sd.beta","Tau2","Alpha","Sd.alpha","Lambda2")
    }
    if(Methodcode==2){
      Result[8]<-NULL
      names(Result)<-c("LB","ResidualVar","Beta","Sd.beta","Tau2","Alpha","Sd.alpha","Delta2","Eta2")
    }
    if(Methodcode==3){
      Result[9]<-NULL; Result[9]<-NULL;
      names(Result)<-c("LB","ResidualVar","Beta","Sd.beta","Sigma2","Alpha","Sd.alpha","Gamma")
    }
    if(Methodcode==4|Methodcode==5){
      Result[9]<-NULL; Result[9]<-NULL;
      names(Result)<-c("LB","ResidualVar","Beta","Sd.beta","Sigma2","Alpha","Sd.alpha","Rho")
      Result$Sigma2<-Result$Sigma2[1]
    }
    if(Methodcode==6){
      Result[9]<-NULL; Result[9]<-NULL;
      names(Result)<-c("LB","ResidualVar","Beta","Sd.beta","Sigma2","Alpha","Sd.alpha","Rho")
      Result$Sigma2<-Result$Sigma2[1:2]
    }
    if(Methodcode==7){
      Result[9]<-NULL; Result[9]<-NULL;
      names(Result)<-c("LB","ResidualVar","Beta","Sd.beta","Sigma2","Alpha","Sd.alpha","Rho")
    }  
    Result$LB<-Result$LB[which(Result$LB!=0)]
    Result$ResidualVar<-Result$ResidualVar[Result$ResidualVar>0]   
    Result
  }
  
  tuning<-function(Pheno,Geno,Methodcode,HyperparametersMatrix,Covariates,Threshold,Maxiterations,Randomize,CVFoldTuning,Printinfo){
    
    N<-length(Pheno)
    Nset<-nrow(HyperparametersMatrix)
    MSE<-numeric(Nset)
    Am<-N%%CVFoldTuning
    if(Am==0){N2<-N}else{N2<-N+CVFoldTuning-Am}
    Shuffled<-matrix(sample(1:N2,N2,replace=FALSE),ncol=CVFoldTuning)
    Test<-as.list(numeric(CVFoldTuning))
    for(fold in 1:CVFoldTuning){Test[[fold]]<-Shuffled[,fold][Shuffled[,fold]<=N]}  
    for(set in 1:Nset){
      Prediction<-numeric(N)
      for(fold in 1:CVFoldTuning){
        if(Printinfo) cat("tuning","set",set,"fold",fold,"\n")
        Result<-genomewideregression(Pheno[-Test[[fold]]],Geno[-Test[[fold]],,drop=FALSE],Methodcode,HyperparametersMatrix[set,],
                                     Covariates[-Test[[fold]],,drop=FALSE],Threshold,Maxiterations,Randomize)
        Prediction[Test[[fold]]]<-Geno[Test[[fold]],,drop=FALSE]%*%matrix(Result$Beta,ncol=1)+Covariates[Test[[fold]],,drop=FALSE]%*%matrix(Result$Alpha,ncol=1)
      }
      MSE[set]<-sqrt(sum((Prediction-Pheno)^2,na.rm=T)/N)
    }
    MSE
  }
  
  #Pheno
  stopifnot(is.vector(Pheno))
  Pheno[is.na(Pheno)]<-99999999 #missing value in the C script
  N<-length(Pheno)
  
  #Geno
  stopifnot(is.matrix(Geno))
  if(any(is.na(Geno))) stop("NA in Geno is not allowed")
  stopifnot(nrow(Geno)==length(Pheno))
  P<-ncol(Geno)
  
  #Method
  Method<-match.arg(Method)
  if(Method=="BL"){
    Methodcode<-1
  }else{
    if(Method=="EBL"){
      Methodcode<-2
    }else{
      if(Method=="wBSR") {
        Methodcode<-3
      }else{
        if(Method=="BayesC"){
          Methodcode<-4
        }else{
          if(Method=="SSVS"){
            Methodcode<-5
          }else{
            if(Method=="MIX"){
              Methodcode<-6
            }else{
              if(Method=="BayesB"){
                Methodcode<-7
              }else{
                stop("Method specification error")
            }}}}}}}
  if(Methodcode==1|Methodcode==2){Priortype<-1} else {Priortype<-2}
  
  #Hyperparameters
  Nh<-c(2,4,3,3,4,4,3)#number of hyperparameters for each method
  stopifnot(is.matrix(Hyperparameters)|is.vector(Hyperparameters))
  if(is.matrix(Hyperparameters)){
    if(ncol(Hyperparameters)!=Nh[Methodcode]){
      stop("Check the number of columns of Hyperperameters")
    }  
    if(Methodcode==1|Methodcode==2)
      if(any(Hyperparameters<=0))stop("Hyperparameters should be positive")
    if(Methodcode==3|Methodcode==4|Methodcode==7){
      if(any(Hyperparameters[,1]<=0))stop("Nu should be >0")
      if(any(Hyperparameters[,2]<0))stop("S2 should be >=0")
      if(any(Hyperparameters[,3]>1|Hyperparameters[,3]<=0))stop("Kappa should be 0<Kappa<=1")
    }
    if(Methodcode==5|Methodcode==6){
      if(any(Hyperparameters[,1]<=0))stop("c should be >0")
      if(any(Hyperparameters[,2]<=0))stop("Nu should be >0")
      if(any(Hyperparameters[,3]<0))stop("S2 should be >=0")  
      if(any(Hyperparameters[,4]>1|Hyperparameters[,4]<=0))stop("Kappa should be 0<Kappa<=1")
    }
    Nset<-nrow(Hyperparameters)
    if(Nset==1)Hyperparameters<-as.vector(Hyperparameters)
  } else {
    if(length(Hyperparameters)!=Nh[Methodcode]){
      stop("Check the length of Hyperperameters")
    }
    if(Methodcode==1|Methodcode==2)
      if(any(Hyperparameters<=0))stop("Hyperparameters should be positive")
    if(Methodcode==3|Methodcode==4|Methodcode==7){
      if(Hyperparameters[1]<=0)stop("v should be >0")
      if(Hyperparameters[2]<0)stop("S2 should be >=0")
      if(Hyperparameters[3]>1|Hyperparameters[3]<=0) stop("Kappa should be 0<Kappa<=1")
    }
    if(Methodcode==5|Methodcode==6){
      if(Hyperparameters[1]<=0)stop("c should be >0")
      if(Hyperparameters[2]<=0)stop("v should be >0")
      if(Hyperparameters[3]<0)stop("S2 should be >=0")  
      if(Hyperparameters[4]>1|Hyperparameters[4]<=0) stop("Kappa should be 0<Kappa<=1")
    }
    Nset<-1
  }

  #Function
  if(Function=="cv"){
    if(Nfold>1){     
      Nfold<-round(Nfold)
      Am<-N%%Nfold
      if(Am==0){N2<-N}else{N2<-N+Nfold-Am}
      Partition<-matrix(sample(1:N2,N2,replace=FALSE),ncol=Nfold)
      Partition[Partition>N]<- -9
      Random<-TRUE        
    }else{
      if(Nfold==-1){
        Partition<-matrix(1:N,ncol=N)
        Nfold<-N
        Random<-FALSE       
      }else{
        if(Nfold==-9){
          if(is.null(Partition)){
            stop("Partition should be specified when Nfold==-9")          
          }else{
            if(!is.matrix(Partition)|!is.numeric(Partition)|any(Partition>N,na.rm=TRUE)|any(is.na(Partition))|any(Partition==0,na.rm=TRUE)){
              stop("Partition matrix error")
            }else{
              Nfold<-ncol(Partition)
              Random<-FALSE
            }                        
          }         
        }else{
          stop ("Nfold specification error")
        }         
      }
    }
  }else{
    if(Function!="fitting"&Function!="tuning")
      stop("Function specification error")
  }
      
  #Covariates
  if(is.character(Covariates)){
    if(Covariates=="Intercept"){
      Covariates<-matrix(rep(1,N),ncol=1); Fc<-1
    } else {
      stop("Covariate specification error")
    }     
  } else {
    stopifnot(is.matrix(Covariates))
    if(any(is.na(Covariates))) stop("NA in Covariates is not allowed")
    stopifnot(nrow(Covariates)==length(Pheno))
    Fc<-ncol(Covariates)
  }
  
  #Max iterations
  stopifnot(Maxiterations>0)
  
  #RandomIni
  if(RandomIni) {Randomize=1} else {Randomize=0}
 
  #Print information
  if(Printinfo){
    cat("\n")
    cat("Method:", Method,"\n")
    cat("Hyperparameters:\n")
    if(Methodcode==1){cat(" Phi","Omega\n")}
    if(Methodcode==2){cat(" Phi","Omega","Psi","Theta\n")} 
    if(Methodcode==3|Methodcode==4|Methodcode==7){cat(" v","S2","Kappa\n")}
    if(Methodcode==5|Methodcode==6){cat(" c","v","S2","Kappa\n")}
    if(Nset==1){
      cat("",Hyperparameters,"\n")      
    }else{
      if(Function=="fitting"){
        cat("",Hyperparameters[1,],"\n")
        cat(" Only the first set is used\n")
      }else{
        for(set in 1:Nset){cat("",Hyperparameters[set,],"\n")}        
      }
    }
    cat("N. of individuals:",N,"\n")
    cat("N. of markers:",P,"\n")
    cat("N. of covariates:",Fc,"\n")
    cat("\n")
  }

  #Calculation
  if(Function=="fitting"){
    if(Printinfo){cat("Model fitting\n")}  
    Result<-genomewideregression(Pheno,Geno,Methodcode,Hyperparameters,Covariates,Threshold,Maxiterations,Randomize)
  }

  if(Function=="tuning"){
    if(Printinfo){cat("Model fitting after hyperparameter tuning\n")}
    if(Nset==1){
      Result<-genomewideregression(Pheno,Geno,Methodcode,Hyperparameters,Covariates,Threshold,Maxiterations,Randomize)      
    }else{
      MSE<-tuning(Pheno,Geno,Methodcode,Hyperparameters,Covariates,Threshold,Maxiterations,Randomize,CVFoldTuning,Printinfo)
      bestset<-which.min(MSE)
      Result<-genomewideregression(Pheno,Geno,Methodcode,Hyperparameters[bestset,],Covariates,Threshold,Maxiterations,Randomize)
    }
    if(Nset>1){
      Result$MSE<-cbind(1:Nset,Hyperparameters,MSE)
      if(Methodcode==1){Hypnames<-c("Phi","Omega")}
      if(Methodcode==2){Hypnames<-c("Phi","Omega","Psi","Theta")} 
      if(Methodcode==3|Methodcode==4|Methodcode==7){Hypnames<-c("Nu","S2","Kappa")}
      if(Methodcode==5|Methodcode==6){Hypnames<-c("c","Nu","S2","Kappa")}    
      colnames(Result$MSE)<-c("Set",Hypnames,"MSE")
      Result$MSE<-data.frame(Result$MSE)
    }         
  }
    
  if(Function=="cv"){
    if(Printinfo){cat("Cross-validation\n")}      
    Prediction<-NULL
    MSE<-NULL
    for(fold in 1:Nfold){
      if(Printinfo){cat("CV fold",fold,"\n")}
      Test<-Partition[,fold]
      Test<-Test[Test!=-9]
      if(Nset==1){
        Result<-genomewideregression(Pheno[-Test],Geno[-Test,,drop=FALSE],Methodcode,Hyperparameters,Covariates[-Test,,drop=FALSE],Threshold,Maxiterations,Randomize)         
      }else{
        MSE.fold<-tuning(Pheno[-Test],Geno[-Test,,drop=FALSE],Methodcode,Hyperparameters,Covariates[-Test,,drop=FALSE],Threshold,Maxiterations,Randomize,CVFoldTuning,Printinfo)
        bestset<-which.min(MSE.fold)
        Result<-genomewideregression(Pheno[-Test],Geno[-Test,,drop=FALSE],Methodcode,Hyperparameters[bestset,],Covariates[-Test,,drop=FALSE],Threshold,Maxiterations,Randomize)        
      }

      BV<-Geno[Test,,drop=FALSE]%*%matrix(Result$Beta,ncol=1)
      Prediction<-rbind(Prediction,
                        cbind(c(1:N)[Test],
                              Pheno[Test],
                              BV+Covariates[Test,,drop=FALSE]%*%matrix(Result$Alpha,ncol=1),
                              BV
                        )
      )
      if(Nset>1){
        MSE<-rbind(MSE,c(fold,bestset,Hyperparameters[bestset,],MSE.fold[bestset]))
      }
    }
    colnames(Prediction)<-c("Test","Y","Yhat","BV")
    Prediction<-data.frame(Prediction)
    if(Nset>1){
      if(Methodcode==1){Hypnames<-c("Phi","Omega")}
      if(Methodcode==2){Hypnames<-c("Phi","Omega","Psi","Theta")} 
      if(Methodcode==3|Methodcode==4|Methodcode==7){Hypnames<-c("Nu","S2","Kappa")}
      if(Methodcode==5|Methodcode==6){Hypnames<-c("c","Nu","S2","Kappa")}
      colnames(MSE)<-c("Fold","ChosenSet",Hypnames,"MSE")
      MSE<-data.frame(MSE)
    }
    if(Random) {Result<-list(Prediction=Prediction,MSE=MSE,Partition=Partition)}else{Result<-list(Prediction=Prediction,MSE=MSE)}    
  }

  if(Printinfo) cat("Finished\n")
  Result
}
