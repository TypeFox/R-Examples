stepMLRMPA <-
function(tree,Clusters,N,op1,op2,tr.tst){
  xtr<-as.data.frame(tr.tst[[4]])
  ytr<-as.data.frame(tr.tst[[1]])
  xtst<-as.data.frame(tr.tst[[5]])
  ytst<-as.data.frame(tr.tst[[2]])
  CutreeVar <- cutreevar(tree,Clusters,matsim=FALSE)
  Cluster <- as.data.frame(CutreeVar$cluster)
  Clusterth_size <- CutreeVar$size
  a<-vector()
  for (i in 1:Clusters){
    a[i]<-paste("Cluster",i,"th",":",Clusterth_size[i],";")
  }
  a2<-t(as.matrix(a))
  a3<-paste("the","result","of","stepwise-MLRMPA")
  a4<-paste("the dim of var selection",":",tr.tst[[5]][1])
  a5<-paste("the dim of cor selection",":",tr.tst[[5]][2])
  write.table(a3,file=op2,append=TRUE,col.names=F,row.names=F,quote=F)
  write.table(a2,file=op2,append=TRUE,col.names=F,row.names=F,quote=F)
  write.table(a4,file=op2,append=TRUE,col.names=F,row.names=F,quote=F)
  write.table(a5,file=op2,append=TRUE,col.names=F,row.names=F,quote=F)

  colname <-t(c("NO.Var","R2.tr","Radj2.tr","LooCV","LgoCV","S2","F_statistic","R2.tst","Q2.ext"))
  write.table(colname,file=op1,append=TRUE,col.names=F,row.names=F)
  Clusterth <- list()
  Var <- list()
  g<-vector()
  xtr<-as.data.frame(xtr)
  for (i in 1:N){
    for (j in 1:Clusters){
      Clusterth[[j]]<-xtr[which(Cluster==j)]
      Var[[j]]<-Clusterth[[j]][sample(dim(Clusterth[[j]])[2],1,replace=T)]
      g[j]<-colnames(Var[[j]])
    }
    for (j in 1:Clusters){
	  if (j==1){
	    matrix<-Var[[1]]
	  }
	  else{
	    matrix<-cbind(matrix,Var[[j]])
	  }
	}
    Var_Y<-cbind(matrix,ytr)
    lm.model <- lm(expr~.,data=Var_Y)
    step.model <- step(lm.model,direction="both",trace=0)
    if (dim(step.model$model)[2]!=1){
      tr.r2   <-summary(step.model)$r.squared
      tr.radj2<-summary(step.model)$adj.r.squared
      LooCV<-lmcv(step.model)[[1]]
      LgoCV<-lmcv(step.model,ng=5)[[1]]
      S2<-summary(step.model)[[6]]^2
      F_statistic <- summary(step.model)$fstatistic[1]    
    
      y.pred.lm<-predict(lm.model,newdata=xtst)
      y2.lm <-y.pred.lm                                                  
      R2.tst.lm <-Rsquared(y2.lm,ytst[,1])
    
      y.pred.slm<-predict(step.model,newdata=xtst)
      y2.slm <-y.pred.slm                                                
      R2.tst <-Rsquared(y2.slm,ytst[,1])
      Q2.ext<-Qsquared.ext(y.pred.slm,ytst[,1],ytr[,1])
      NO.Var<-length(which(strsplit(as.character(summary(step.model)$call)[2]," ")[[1]]=="+"))+1
      
      result<-t(as.matrix(c(NO.Var,tr.r2,tr.radj2,LooCV,LgoCV,S2,F_statistic,R2.tst,Q2.ext)))
      write.table(result,file=op1,append=TRUE,col.names=F,row.names=F)       
    
      step.model1<-matrix(0,6,1)
      step.model1[1,1]<-paste("model:clusters",Clusters,"          the order of samples " ,i)
      step.model1[2,1]<-paste("multiple","linear","regression","model")
      step.model1[3,1]<-paste("y","~",paste(rownames(summary(lm.model)$coef)[1:Clusters+1],collapse="+"))
      step.model1[4,1]<-paste("tr_r.squared",":",summary(lm.model)$r.squared)
      step.model1[5,1]<-paste("tr_adj.r.squared",":",summary(lm.model)$adj.r.squared)
      step.model1[6,1]<-paste("tst_r.squared",":",R2.tst.lm)
      write.table(step.model1,file=op2,append=TRUE,col.names=F,row.names=F,quote=F)    
    
      colnames.lm<-matrix(0,1,5)
      colnames.lm[1,1]<-c("variables")
      colnames.lm[1,2:5]<-colnames(summary(step.model)$coefficients)
      write.table(colnames.lm,file=op2,append=TRUE,col.names=F,row.names=F,quote=F)
      write.table(summary(lm.model)$coefficients,file=op2,append=TRUE,col.names=F,quote=F)
      write.table("\n",file=op2,append=TRUE,col.names=F,row.names=F,quote=F)
      step.model2<-matrix(0,8,1)
      step.model2[1,1]<-paste("multiple","linear","stepwise","regression","model")
      step.model2[2,1]<-as.character(summary(step.model)$call)[2]
      step.model2[3,1]<-paste("tr_r.squared",":",summary(step.model)$r.squared)
      step.model2[4,1]<-paste("tr_adj.r.squared",":",summary(step.model)$adj.r.squared)
      step.model2[5,1]<-paste("tst_.r.squared",":",R2.tst)
      step.model2[6,1]<-paste("LooCV",":",LooCV)
      step.model2[7,1]<-paste("LgoCV",":",LgoCV)
      step.model2[8,1]<-paste("tst_r.squared.cv",":",Q2.ext)
      write.table(step.model2,file=op2,append=TRUE,col.names=F,row.names=F,quote=F)
      colnames<-matrix(0,1,5)
      colnames[1,1]<-c("variables")
      colnames[1,2:5]<-colnames(summary(step.model)$coefficients)
      write.table(colnames,file=op2,append=TRUE,col.names=F,row.names=F,quote=F)
      write.table(summary(step.model)$coefficients,file=op2,append=TRUE,col.names=F,quote=F)
      write.table("\n",file=op2,append=TRUE,col.names=F,row.names=F,quote=F)      
	  }
	  if(i==1){
          print(paste(Clusters,":" ,i,Sys.time()))
          time1<-Sys.time()
      	  }
	    if(i%%100==0 & i!=N){
          print(paste(Clusters,":" ,i,Sys.time()))
     	  }
         if(i==N){
          print(paste(Clusters,":" ,i,Sys.time()))
          time2<-Sys.time()
          print(difftime(time2,time1,units=c("mins")))
	     }
   }
  return
  Clusterth_size  
  }
