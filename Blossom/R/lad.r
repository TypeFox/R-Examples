lad<-function(formula,data,contrasts=NULL,number.perms=5000, quant,test=FALSE,all.quants=FALSE,OLS,weights){
      #some default settings for lad regression with no hypothesis testing and such
      i_Test=1
      i_NumToDrop=0 #zero for regression something else for hypothesis testing
      i_LaHy=1
      l_DoRankScore=FALSE
      l_HasHypCmdLine=FALSE
      l_IsLaPg=TRUE
      if(missing(OLS)) OLS=FALSE
      if(OLS) l_IsLaPg=FALSE
      l_IsLaPgType=TRUE
      i_LaPgType=1
      if(missing(quant)) theta=-.5
      else theta=quant

 ##################################
       Call<-match.call()
      if (missing(data))
              data <- environment(formula)
      mf <- match.call(expand.dots = FALSE)
      m <- match(c("formula", "data","weights"), names(mf), 0L)
       mf <- mf[c(1L, m)]
          mf$drop.unused.levels <- TRUE
          mf[[1L]] <- as.name("model.frame")
          mf <- eval(mf, parent.frame())          
       weights <- as.vector(model.weights(mf))
      y <- model.response(mf, "any")
          if (length(dim(y)) == 1L) {
              nm <- rownames(y)
              dim(y) <- NULL
              if (!is.null(nm))
                  names(y) <- nm
          }
          mt <- attr(mf, "terms")
          
          x <- if (!is.empty.model(mt))
                     model.matrix(mt, mf, contrasts)
           contr <- attr(x, "contrasts")     
         if(!is.null(weights)){
             y<-y*weights
             x<-x*weights
           } 
      #Number of variables doesn't include an intercept so remove it if present
      NumVars<-dim(x)[2]-attr(mt,"intercept")
      la_HasIntercept=rep(as.logical(attr(mt,"intercept")),times=2)
      NumObs=dim(x)[1]

      input.x<-x
      #remove incomplete cases
      i_NumObs=nrow(x)
         if(ncol(x)==1) x<-na.omit(x)
            else x<-x[complete.cases(x),]
        comp.cases<-nrow(x)
        if(comp.cases/i_NumObs!=1) {warning(paste(i_NumObs-comp.cases,ifelse(i_NumObs-comp.cases==1," case was removed because of a missing value ",
          " cases were removed because of missing values."),sep=""))
             NumObs<-nrow(x)   
          }
        
      response<-matrix(data=y,ncol=1)
      colnames(response)<-names(mf[1])

      if(!all.quants & !OLS){
        x1<-matrix(data=0,nrow=(dim(x)[1]+2),ncol=dim(x)[2]+2)
        x1[1:dim(x)[1],1:dim(x)[2]]<-x
        x1[1:dim(x)[1],(dim(x)[2]+1)]<-y
        x1[,dim(x1)[2]]<-seq(1:dim(x1)[1])
        x<-x1

         i_NumRows=NumObs+2
      } else {
        x<-cbind(x,y)
        i_NumRows=NumObs
        }

      
      #think about how to rearange the data later
      #x1<-matrix(data=NA,nrow=(NumObs+2),ncol=(NumVars+3)) #for regression we get two extra rows and columns plus a contant column for intercept
      if(all.quants) i_NumRows=dim(x)[1]
      i_NumCols=dim(x)[2]
      
      da_Betas<-rep(0,times=(NumVars+attr(mt,"intercept")))
      da_RedBetas<-rep(0,times=(NumVars+attr(mt,"intercept")-i_NumToDrop))

      #OLS counts the intercept as a variable
      if(OLS) NumVars<-NumVars+attr(mt,"intercept")
      
       ia_PosToDrop<-rep(0,times=(NumVars+la_HasIntercept[1]))
      da_STV<-0


      i_NumRedVars<-NumVars-i_NumToDrop
      ia_NumLaVars<-c(NumVars,i_NumRedVars)
      if(OLS) ia_NumLaVars<-ia_NumLaVars-la_HasIntercept
      da_ResRed<-rep(0,times=NumObs)
      da_Resids<-rep(0,times=NumObs)
      da_Sol<-matrix(data=0,ncol=(4*NumObs*all.quants),nrow=((NumVars+3)*all.quants+1))

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
       i_NumVarRed<-ia_NumLaVars[1]+la_HasIntercept[1]
      d_To<-d_Tn_RS<-i_CntSumAVR<-d_PValue<-d_PVTn<-d_SumAbsValRes<-d_SumAbsValResRed<-d_WtSumAbsDevsFulMod<-d_WtSumAbsDevsRedMod<-i_Iter<-i_ExitCode<-iEr<-i_lsol<-0
      ladOut=.Fortran("wraplad",
                              x,
                              as.integer(NumObs),
                              as.integer(NumVars),
                              as.double(theta),
                              as.integer(number.perms),
      	                      as.integer(i_NumToDrop),
                              ia_PosToDrop,
                              as.numeric(0),
                              as.integer(i_Test),
      	                      as.logical(FALSE),
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
                              as.integer(all.quants),
                              as.logical(test),
                              as.logical(l_DoRankScore),
                              as.logical(OLS),
                              as.logical(l_HasHypCmdLine),
                              as.logical(l_IsLaPg),
                              as.integer(i_LaPgType),
                              ia_NumLaVars,
                              da_ResRed,
                              da_Resids,
                              as.integer(i_NumRows),
                              as.integer(i_NumCols),
                              as.integer(ia_NumLaVars[1]+la_HasIntercept[1]),
                              as.integer(i_NumVarRed),
                              da_Sol,
                              as.integer(i_lsol))

                                  QuantOut=matrix(0,0,0)
                                 if(all.quants==TRUE) QuantOut=ladOut[[42]][1:(NumVars+4),1:ladOut[[43]]]
          #here I set reduced output equal to full output to streamline since here the
          #reduced and the full model are the same and it streamlines hypothesis testing
        LadOut<-new("LADObj",inputData=input.x,NumObs=ladOut[[2]],NumVars=ladOut[[3]],
          theta=ladOut[[4]], NumPerm=ladOut[[5]],
          Test=ladOut[[9]],DoublePermutation=ladOut[[10]],T_o=ladOut[[11]],
          AsyRankScore=ladOut[[12]],P_value=ladOut[[14]],P_valueTN=ladOut[[15]],
          Betas=ladOut[[16]],RedBetas=ladOut[[16]],SumAbsValRes=ladOut[[18]], SumAbsValResRed=ladOut[[18]],
          WtSumAbsDevsFulMod=ladOut[[20]], WtSumAbsDevsRedMod=ladOut[[21]],
      	  NumIter=ladOut[[22]],ExitCode=ladOut[[23]],PermVals=ladOut[[25]],
          HasIntercept=ladOut[[27]],DoAllQuants=ladOut[[28]],
          DoRankScore=ladOut[[30]],IsOLS=ladOut[[31]],
          NumLaVars=ladOut[[35]],ResRed=ladOut[[36]],Resids=ladOut[[37]],
          Call=Call,response=response,QuantOut=QuantOut)

LadOut
}
