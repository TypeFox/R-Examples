tri.app <-
function(ms,ET,M.exp,E.exp,T.exp,N = 0.25,method="pearson",iqr.filter = c(log2(1.5),log2(1.5),log2(1.5)),
                  cor.MvsET = c(0.3,0.3),cor.EvsT.dif = 0.45,cor.EvsT.abs = 0.4,
                  ET.fc.filter = log2(1.5),ET.p.filter = 0.01,rand=100,correction="BH",cores=1){
  #preprocess the mRNA and lncRNA expression profile by IQR and filter the TF-target regulations in which the TF(or target) expression is not in the preprocessed expression profile
  iqr <- apply(M.exp, 1, IQR, na.rm=TRUE) 
  data.M <- M.exp[iqr>iqr.filter[1],] 
  iqr <- apply(E.exp, 1, IQR, na.rm=TRUE) 
  data.E <- E.exp[iqr>iqr.filter[2],] 
  iqr <- apply(T.exp, 1, IQR, na.rm=TRUE) 
  data.T <- T.exp[iqr>iqr.filter[3],]  
  index.ET.E<-as.character(ET[,1])%in%rownames(data.E)
  index.ET.T<-as.character(ET[,2])%in%rownames(data.T)
  index.ET<-index.ET.E&index.ET.T
  ET<-ET[index.ET,]
 
  if(dim(ET)[1]==0){
  tmptmpMET<-list()
  tmptriplets<-as.data.frame(matrix(nrow=0,ncol=7))
  colnames(tmptriplets)<-c("modulator","effector","target","R_low","R_HIGH","p-value","fdr")
  tmptmpMET[["triplets"]]<-tmptriplets
  tmpfactor<-"aa"
  tmpfactorindex<-F
  tmptmpMET[["initialnot"]]<-tmpfactor[tmpfactorindex]
  tmptmpMET[["filterdnot"]]<-tmpfactor[tmpfactorindex]
  return(tmptmpMET)
  }else{
  #predict triplet by lncRNAs(one by one) usiong the mclapply function in parallel package to improve speed
  #define a internal function for predicting
  infun<-function(m){
    if(m%in%rownames(data.M)){
    #an object named "MET" was constructed to contain the predicted triplets
    MET<-as.data.frame(matrix(nrow=0,ncol=7))
    colnames(MET)<-c("modulator","effector","target","low_correlation","high_correlation","p-value","fdr")

file <- paste("correlation_matrix_",m,".Rdata",sep="")
      if(file.exists(file)){
        load(file)
      }else{
        
        exp.threshold <- quantile(data.M[m,],probs=N)[1,1]
        LOW.index <- data.M[m,]<exp.threshold
        LOW<- apply(ET,1,function(x) cor(t(data.T[as.character(x[2]),LOW.index]),t(data.E[as.character(x[1]),LOW.index]),use="pairwise.complete.obs",method=method))
        exp.threshold <- quantile(data.M[m,],probs=1-N)[1,1]
        HIGH.index <- data.M[m,]>exp.threshold
        HIGH<- apply(ET,1,function(x) cor(t(data.T[as.character(x[2]),HIGH.index]),t(data.E[as.character(x[1]),HIGH.index]),use="pairwise.complete.obs",method=method))
      }
      
      data.EE<-data.E[as.character(ET[,1]),]
      DE <- apply(data.EE, 1, function(x) t.test(x[LOW.index], x[HIGH.index])) 
      DE.p <- sapply(DE, "[[", "p.value")  
      DE.FC <- sapply(DE, function(x) diff(x$estimate)) 
      DE.FC<-abs(DE.FC)
      index.DE <- !(DE.p < ET.p.filter & DE.FC > ET.fc.filter) 
      
      data.MvsT.cor <-cor(t(data.T[as.character(ET[,2]),]),t(data.M[m,]))
      data.MvsE.cor <- cor(t(data.E[as.character(ET[,1]),]),t(data.M[m,]))
      index.MvsT<-data.MvsT.cor<cor.MvsET[2]
      index.MvsE<-data.MvsE.cor<cor.MvsET[1]
      index.MvsET<-index.MvsT&index.MvsE
      
      delta.cor<-HIGH-LOW
      index.diff <- abs(delta.cor)>cor.EvsT.dif 
      index.abs <- abs(LOW)>cor.EvsT.abs | abs(HIGH)>cor.EvsT.abs 
      index <- index.diff & index.abs 
      
      index<-index&index.MvsET&index.DE
      if(sum(index)!=0){
        tmp<-cbind(rep(m,sum(index)),ET[index,],LOW[index],HIGH[index])
        
        ml<-dim(data.M)[2]
        LOW.indexs<-sum(LOW.index)
        HIGH.indexs<-sum(HIGH.index)
        deltar<-tmp[5]-tmp[4]
        for(rr in 1:rand){
          LOW.rand.index<-sample(1:ml,LOW.indexs)
          HIGH.rand.index<-sample(setdiff(1:ml,LOW.rand.index),HIGH.indexs)
          LOW.rand<- apply(tmp[,2:3],1,function(x) cor(t(data.T[as.character(x[2]),LOW.rand.index]),t(data.E[as.character(x[1]),LOW.rand.index]),use="pairwise.complete.obs",method=method))
          HIGH.rand<- apply(tmp[,2:3],1,function(x) cor(t(data.T[as.character(x[2]),HIGH.rand.index]),t(data.E[as.character(x[1]),HIGH.rand.index]),use="pairwise.complete.obs",method=method))
          deltar<-cbind(deltar,HIGH.rand-LOW.rand)
        }
        p<-numeric()
        for(r in 1:dim(deltar)[1]){
          if(deltar[r,1]<=0){tmpp<-sum(deltar[r,2:(rand+1)]<deltar[r,1])/rand}
          if(deltar[r,1]>0){tmpp<-sum(deltar[r,2:(rand+1)]>deltar[r,1])/rand}
          p<-c(p,tmpp)
        }
        fdr<-p.adjust(p,method=correction)
        tmp<-cbind(tmp,p,fdr)
        colnames(tmp)<-c("modulator","effector","target","R_low","R_HIGH","P_value","fdr")
MET<-rbind(MET,tmp)
}
    }else{
      if(m%in%rownames(M.exp)){
    MET<-"filterdnot"
      }else{
  MET<-"initialnot"
      }
    }
    return(MET)
  }#for infun
  
  #producing result
  tmpMET<-mclapply(ms,infun,mc.cores=cores)
  names(tmpMET)<-ms
  #normalized format
  findex<-grep("filterdnot",tmpMET)
  iindex<-grep("initialnot",tmpMET)
  index<-setdiff(1:length(tmpMET),c(findex,iindex))
  tmptmpMET<-list()
  tmptri<-as.data.frame(matrix(nrow=0,ncol=7))
  colnames(tmptri)<-c("modulator","effector","target","R_low","R_HIGH","p-value","fdr")
  for(tmpi in index){
  tmptri<-rbind(tmptri,tmpMET[[tmpi]])
  }
  tmptmpMET[["triplets"]]<-tmptri
  tmptmpMET[["initialnot"]]<-ms[iindex]
  tmptmpMET[["filterdnot"]]<-ms[findex]
  
  return(tmptmpMET)
  }
}
