hypothesis.test<-function(object1,object2,number.perms=5000,save.test=TRUE,double.permutation=FALSE,rank.score=FALSE){
      #some default settings for lad regression with no hypothesis testing and such

       Call<-match.call()

       i_Test=2
      i_NumToDrop=0 #zero for regression something else for hypothesis testing
      i_LaHy=2
      l_DoAllQuants=FALSE
      if(object1@IsOLS!=object2@IsOLS) stop("Models must be fit using the same method (OLS or LAD)")
      l_IsOLS=object1@IsOLS
      l_HasHypCmdLine=TRUE

      if(l_IsOLS){
        l_IsLaPg=TRUE
        l_IsLaPgType=FALSE
        i_LaPgType=2
        test=FALSE}
      else{
        l_IsLaPg=FALSE
        l_IsLaPgType=TRUE
        i_LaPgType=1
        test=TRUE}
  #assume object1 is the full model  and if object2 is the full model switch them

 if(sum(is.na(match(colnames(object2@inputData),colnames(object1@inputData))))>0){
  temp<-object1
  object1<-object2
  object2<-temp
  #get columns in the same order
  #column.names<-colnames(object1@inputData)[(match(colnames(object2@inputData),colnames(object1@inputData)))]
  #object1@inputData<-as.matrix(object1@inputData[,match(colnames(object2@inputData),colnames(object1@inputData))])
  }

  #now all of object2 names should be in object1 otherwise the models aren't nested
  #check all names are the same and that the corresponding data is the same
    if(any(object1@inputData[,match(colnames(object2@inputData),
    colnames(object1@inputData))]!=object2@inputData) | any(object1@response!=object2@response)) stop("Models are not nested")

    if(object1@theta!=object2@theta) stop("Different quantiles were used for the two model fits")
      else theta<-object1@theta

 ia_PosToDrop<-as.numeric(is.na(match(colnames(object1@inputData),colnames(object2@inputData))))
 i_NumToDrop<-sum(ia_PosToDrop)

      #Number of variables doesn't include an intercept so remove it if present
      NumVars<-object1@NumVars
      la_HasIntercept=c(object1@HasIntercept[2],object2@HasIntercept[2])
      NumObs=object1@NumObs
      if(object1@IsOLS) ia_NumLaVars<-c(object1@NumVars,object2@NumVars)-la_HasIntercept

      if(!object1@IsOLS){
        x1<-matrix(data=0,nrow=(dim(object1@inputData)[1]+2),ncol=(dim(object1@inputData))[2]+2)
      x1[1:dim(object1@inputData)[1],1:dim(object1@inputData)[2]]<-object1@inputData
      x1[1:dim(object1@inputData)[1],(dim(object1@inputData)[2]+1)]<-object1@response
      x1[,dim(x1)[2]]<-seq(1:dim(x1)[1])
      x<-x1
      
              
         i_NumRows=object1@NumObs+2
      } else {
        x<-cbind(object1@inputData,object1@response)
        i_NumRows=object1@NumObs
        }
      #think about how to rearange the data later
      #x1<-matrix(data=NA,nrow=(NumObs+2),ncol=(NumVars+3)) #for regression we get two extra rows and columns plus a contant column for intercept

      i_NumCols=dim(x)[2]
      da_Betas<-rep(0,times=length(object1@Betas))
      da_RedBetas<-rep(0,times=length(object1@Betas))

      da_STV<-rep(0,times=number.perms)

    
      i_NumRedVars<-NumVars-i_NumToDrop
      ia_NumLaVars<-c(NumVars,i_NumRedVars)
      if(l_IsOLS) ia_NumLaVars<-rep((NumVars-i_NumToDrop),times=2)
      da_ResRed<-rep(0,times=NumObs)
      da_Resids<-rep(0,times=NumObs)
      da_Sol<-matrix(data=0,ncol=(4*NumObs*l_DoAllQuants),nrow=((NumVars+3)*l_DoAllQuants+1))
      
      storage.mode(x)<-"double"
      storage.mode(da_STV)<-"double"
      storage.mode(ia_PosToDrop)<-"integer"
      storage.mode(da_Betas)<-"double"
      storage.mode(da_RedBetas)<-"double"
      storage.mode(ia_NumLaVars)<-"integer"
      storage.mode(la_HasIntercept)<-"logical"
      storage.mode(da_ResRed)<-"double"
      storage.mode(da_Resids)<-"double"
      storage.mode(da_Sol)<-"double"
      
      d_To<-d_Tn_RS<-i_CntSumAVR<-d_PValue<-d_PVTn<-d_SumAbsValRes<-d_SumAbsValResRed<-d_WtSumAbsDevsFulMod<-d_WtSumAbsDevsRedMod<-i_Iter<-i_ExitCode<-iEr<-i_lsol<-0

      ladOut=.Fortran("wraplad",
                              x,
                              as.integer(NumObs),
                              as.integer(NumVars),
                              as.double(theta),
                              as.integer(number.perms),
      	                      as.integer(i_NumToDrop),
                              ia_PosToDrop,
                              as.integer(save.test),
                              as.integer(i_Test),
      	                      as.logical(double.permutation),
                              as.double(d_To),
                              as.double(d_Tn_RS),
                              as.integer(i_CntSumAVR),
                              as.double(d_PValue),
                              as.double(d_PVTn),
      	                      da_Betas,
                              da_RedBetas,
      	                      as.double(d_SumAbsValRes),
                              as.double(d_SumAbsValResRed),
      	                      as.double(d_WtSumAbsDevsFulMod),
                              as.double(d_WtSumAbsDevsRedMod),
      	                      as.integer(i_Iter),
                              as.integer(i_ExitCode),
                              as.integer(iEr),
                              da_STV,
                              as.integer(i_LaHy),
                              la_HasIntercept,
                              as.integer(l_DoAllQuants),
                              as.logical(test),
                              as.logical(rank.score),
                              as.logical(l_IsOLS),
                              as.logical(l_HasHypCmdLine),
                              as.logical(l_IsLaPg),
                              as.integer(i_LaPgType),
                              ia_NumLaVars,
                              da_ResRed,
                              da_Resids,
                              as.integer(i_NumRows),
                              as.integer(i_NumCols),
                              as.integer(ia_NumLaVars[1]+la_HasIntercept[1]),
                              as.integer(ia_NumLaVars[1]+la_HasIntercept[1]),
                              da_Sol,
                              as.integer(i_lsol))
                            
        LadOut<-new("LADObj",inputData=object2@inputData,NumObs=ladOut[[2]],NumVars=ladOut[[3]],
          theta=ladOut[[4]], NumPerm=ladOut[[5]],
          Test=ladOut[[9]],DoublePermutation=ladOut[[10]],T_o=ladOut[[11]],
          AsyRankScore=ladOut[[12]],P_value=ladOut[[14]],P_valueTN=ladOut[[15]],
          Betas=ladOut[[16]],RedBetas=ladOut[[17]][1:(length(ia_PosToDrop)-sum(ia_PosToDrop))],SumAbsValRes=ladOut[[18]], SumAbsValResRed=ladOut[[19]],
          WtSumAbsDevsFulMod=ladOut[[20]], WtSumAbsDevsRedMod=ladOut[[21]],
      	  NumIter=ladOut[[22]],ExitCode=ladOut[[23]],PermVals=ladOut[[25]],
          HasIntercept=ladOut[[27]],DoAllQuants=ladOut[[28]],
          DoRankScore=ladOut[[30]],IsOLS=ladOut[[31]],
          NumLaVars=ladOut[[35]],ResRed=ladOut[[36]],Resids=ladOut[[37]],
          Call=Call,response=object1@response,
          full.mod.names=colnames(object1@inputData))
          LadOut
    }
