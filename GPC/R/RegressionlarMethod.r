RegressionlarMethod <- function(InputDim,Designs,InputDistrib,ParamDistrib,pmaxi,jmax,Output,PCSpace,
                                   Q2tgt,SeedSob,EpsForw,EpsBack,EnrichStep,p,SaveQmax,SaveAmax)
  {
  #print(c("DesignLength",Designs$DesignLength))
  restartAnalysis<-FALSE    
  #Ap = matrix(rep(0,InputDim),1,InputDim)
  #A = Ap
  
#  A=rbind(A,diag(1,3,3))
#  A=rbind(A,c(1,1,0))
#  A=rbind(A,c(0,3,0))
#  A
  
#  DM=DataMatrix(A,Designs$PCE,InputDistrib,pmaxi,PCSpace)
  
  #R2Q2_0 <- c(Reg$R2,Reg$Q2)
  #Qmax<-0
  #SaveQmax<-c(Qmax,0,0)
  #SaveAmax1=SaveAmax[[1]]
  #SaveAmax2=c()
  #SaveAmax3=c()
  #pmaxi=4
  while( SaveQmax[1]<=Q2tgt & p<pmaxi)
    {
    p=p+1    
    print(c('p is',p))
    ACurrent=HIS(p,InputDim,1)#t(indexCardinal(InputDim,p))
#     print('Nbre de' nrow(Acurrent))
#     print(ACurrent)
    #print(Designs$PCE)
    DM=DataMatrix(ACurrent,Designs$PCE,InputDistrib,pmaxi,PCSpace)
    LARSresult=lars(DM, t(t(Output)), type = c("lar"),trace = FALSE, normalize = TRUE, intercept = TRUE, eps = .Machine$double.eps)
    a=as.numeric(unlist(LARSresult$actions))
    Q=c()
    for(i in 1:length(a))
      {
      Atemp=ACurrent[c(1,a[1:i]),]
      Reg <-RegressionLar(Atemp,Designs$PCE,InputDistrib,pmaxi,Output,PCSpace,p)
      restartAnalysis <- Reg$rAnaly
      if(restartAnalysis==TRUE)
        {
        #print("Regression Enrichment");
        Designs$DesignLength <- Designs$DesignLength + EnrichStep
        ED<-GetDesignReg(Designs$DesignLength,InputDim,SeedSob,InputDistrib,ParamDistrib,PCSpace)
        return(list(Designs=list(PCE=ED$PCE,Physic=ED$Physic,Sobol=ED$Sobol,DesignLength=Designs$DesignLength),
                                 p=0,SaveQmax=c(-100,-100,-100),SaveAmax=list(c(0),c(0),c(0))))
        }
      Q[i]=Reg$Q2
      }
    Qmax=max(Q)
    Amax=ACurrent[c(1,a[1:which.max(Q)]),]
    SaveQmax[3]=SaveQmax[2]
    SaveQmax[2]=SaveQmax[1]
    SaveQmax[1]=Qmax
    SaveAmax[[3]]=SaveAmax[[2]]
    SaveAmax[[2]]=SaveAmax[[1]]
    SaveAmax[[1]]=Amax
    #print(SaveQmax)
    #print(SaveAmax)
    if(SaveQmax[1]<=SaveQmax[2] & SaveQmax[1]<=SaveQmax[3])
      {
      #print("Overfitting Enrichment");
      Designs$DesignLength <- Designs$DesignLength + EnrichStep
      ED<-GetDesignReg(Designs$DesignLength,InputDim,SeedSob,InputDistrib,ParamDistrib,PCSpace)      
      SaveQmax<-c(SaveQmax[3],0,0)
      SaveAmax[[1]]=SaveAmax[[3]]
      SaveAmax[[2]]=c(0)
      SaveAmax[[3]]=c(0)
      return(list(Designs=list(PCE=ED$PCE,Physic=ED$Physic,Sobol=ED$Sobol,DesignLength=Designs$DesignLength),p=(p-2),SaveQmax=SaveQmax,SaveAmax=SaveAmax))
      }    
    } 
  return(list(TruncSet=SaveAmax[[1]],Q2=SaveQmax[1],Designs=list(PCE=Designs$PCE,Physic=Designs$Physic,Sobol=Designs$Sobol,DesignLength=Designs$DesignLength,p=p)))  
  }
#   
#   Reg <- Regression(A,Designs$PCE,InputDistrib,pmaxi,Output,PCSpace)
#   R2Q2_0 <- c(Reg$R2,Reg$Q2)
#   print(p)
#   while( R2Q2_0[2]<=Q2tgt & p<pmaxi){      
#     p = p+1
#     print(p)
#     j = 1
#     App=c()
#     while( R2Q2_0[2]<=Q2tgt & j<=min(p,jmax)){              
#       J_jp=c()
#       TC_jp = t(indexCardinal(InputDim,p)) #HIS(p,InputDim,1)
#       C_jp = matrix(TC_jp[which( rowSums(TC_jp>0)==j & rowSums(TC_jp)==p ),],ncol=InputDim)
#       k=0
#       for(i in 1:nrow(C_jp)){
#         A = rbind(Ap,C_jp[i,])
#         Reg <- Regression(A,Designs$PCE,InputDistrib,pmaxi,Output,PCSpace)
#         R2Q2 <- c(Reg$R2,Reg$Q2)
#         restartAnalysis <- Reg$rAnaly
#         if(restartAnalysis==TRUE){
#           print("Enrich2");
#           Designs$DesignLength <- Designs$DesignLength + EnrichStep
#           ED<-GetDesignReg(Designs$DesignLength,InputDim,SeedSob,InputDistrib,ParamDistrib,PCSpace)
#           return(list(Designs=list(PCE=ED$PCE,Physic=ED$Physic,Sobol=ED$Sobol,DesignLength=Designs$DesignLength)))
#         }
#         DeltaR2 = R2Q2[1]-R2Q2_0[1]
#         if(DeltaR2>=EpsForw){ J_jp = rbind(J_jp,c(C_jp[i,],DeltaR2)); k=k+1 }          
#       }
#       if(k>0 & restartAnalysis==FALSE){
#         J_jpStar=matrix(J_jp[sort(J_jp[,InputDim+1],index.return=TRUE,decreasing=TRUE)$ix,(1:InputDim)],ncol=InputDim)
#         R_jp=c()
#         for(i in 1:nrow(J_jpStar)){        
#           DM=DataMatrix(rbind(Ap,J_jpStar[i,]),Designs$PCE,InputDistrib,pmaxi,PCSpace)
#           IM=t(DM)%*%DM
#           if(rcond(IM)<10^(-4)){
#             print("Enrich1")
#             Designs$DesignLength <- Designs$DesignLength + EnrichStep
#             ED<-GetDesignReg(Designs$DesignLength,InputDim,SeedSob,InputDistrib,ParamDistrib,PCSpace)
#             return(list(Designs=list(PCE=ED$PCE,Physic=ED$Physic,Sobol=ED$Sobol,DesignLength=Designs$DesignLength)))
#           }
#           else { 
#             R_jp = rbind(R_jp,J_jpStar[i,]) 
#           }
#         }
#         App = rbind(Ap,R_jp)
#         Ap=App
#       }
#       j=j+1        
#     }
#     if(restartAnalysis==FALSE){
#       Reg <- Regression(Ap,Designs$PCE,InputDistrib,pmaxi,Output,PCSpace)
#       R2Q2_0 <- c(Reg$R2,Reg$Q2)
#       restartAnalysis <- Reg$rAnaly
#       IndSup=c()
#       for(i in 1:nrow(Ap)){
#         A=matrix(Ap[-i,],ncol=InputDim)
#         Reg <- Regression(A,Designs$PCE,InputDistrib,pmaxi,Output,PCSpace)
#         R2Q2 <- c(Reg$R2,Reg$Q2)
#         restartAnalysis <- Reg$rAnaly
#         if(restartAnalysis==TRUE){
#           print("Enrich3");
#           Designs$DesignLength <- Designs$DesignLength + EnrichStep
#           ED<-GetDesignReg(Designs$DesignLength,InputDim,SeedSob,InputDistrib,ParamDistrib,PCSpace)
#           return(list(Designs=list(PCE=ED$PCE,Physic=ED$Physic,Sobol=ED$Sobol,DesignLength=Designs$DesignLength)))
#         }
#         DeltaR2 = R2Q2_0[1]-R2Q2[1]
#         if(DeltaR2<=EpsBack){ IndSup=cbind(IndSup,i) }
#       }
#       if(length(IndSup)!=0 & restartAnalysis==FALSE){
#         Ap<-Ap[-IndSup,]
#         Reg <- Regression(A,Designs$PCE,InputDistrib,pmaxi,Output,PCSpace)
#         R2Q2_0 <- c(Reg$R2,Reg$Q2)
#         restartAnalysis <- Reg$rAnaly
#       }
#     }      
#   }
#   return(list(TruncSet=Ap,R2=R2Q2_0[1],Q2=R2Q2_0[2],Designs=list(PCE=Designs$PCE,Physic=Designs$Physic,Sobol=Designs$Sobol,DesignLength=Designs$DesignLength)))    
# }  
