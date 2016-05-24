binomialCode <-function(data,alternative="two.sided",npNumbers=100,beta=.001,interval=FALSE,method="Z-pooled",cond.row=TRUE,to.plot=TRUE,ref.pvalue=TRUE){
  
  maxPvalue<-function(tbls){
    #Instead of calculating the binomial probabilities several times for a more extreme cell,
    #calculate the probability once and then combine with other more extreme cells:
    xTbls<-tbls[,1]
    yTbls<-tbls[,2]
    nTbls<-length(xTbls) # Number of 'as or more extreme' tables
    
    x.unique<-unique(xTbls)
    y.unique<-unique(yTbls)
    lxu<-length(x.unique)
    lyu<-length(y.unique)
    
    xnr<-max(xTbls)+1
    A<-matrix(nrow=xnr, ncol=npNumbers)
    cellnr<-rep.int(x.unique+1, npNumbers)+xnr*rep(seq(npNumbers)-1, each=lxu)
    A[cellnr]<-dbinom(rep.int(x.unique, npNumbers), Ns[1], rep(int, each=lxu))
    
    ynr<-max(yTbls)+1
    B<-matrix(nrow=ynr, ncol=npNumbers)
    cellnr<-rep.int(y.unique+1, npNumbers)+ynr*rep(seq(npNumbers)-1, each=lyu)
    B[cellnr]<-dbinom(rep.int(y.unique, npNumbers), Ns[2], rep(int, each=lyu))
    
    prob<-matrix(A[xTbls+1,]*B[yTbls+1,],ncol=length(int))
    prob<-colSums(prob, nTbls)
    
    np<-int[which(prob==max(prob))]
    pvalue<-max(prob)+beta
    
    list<-list(prob=prob,pvalue=pvalue,np=np)
    return(list)
  }
  
  #If conditioning on row, then transpose 2x2 table
  if(cond.row){data<-t(data)}
  
  Ns<-colSums(data,dims=1)
  N<-sum(Ns)
  
  #Specify nuisance parameter range
  if(interval){
    tempInt<-binom.test(sum(data[1,]), N, conf.level=1-beta)$conf.int[1:2]
    int<-seq(max(c(0.00001,tempInt[1])),min(c(0.99999,tempInt[2])),length=npNumbers)
  } else {int<-seq(0.00001,.99999,length=npNumbers); beta<-0}
  
  x <- rep(0:Ns[1], each=(Ns[2]+1))
  y <- rep.int(0:Ns[2], Ns[1]+1)
  p1 <- x/Ns[1]
  p2 <- y/Ns[2]
  
  #Calculate the test statistics:
  if(tolower(method)=="boschloo"){
    TX<-matrix(c(x,y,apply(matrix(c(x,Ns[1]-x,y,Ns[2]-y),(Ns[1]+1)*(Ns[2]+1),4),1,
                           FUN=function(tbls){fisher.2x2(matrix(tbls,2,2),alternative=tolower(alternative))})),
               (Ns[1]+1)*(Ns[2]+1),3)
  }
  
  if(tolower(method) %in% c("z-pooled","pooled","score")){
    TX<-matrix(c(x,y,(p1-p2)/sqrt(((x+y)/N)*(1-((x+y)/N))*sum(1/Ns))),(Ns[1]+1)*(Ns[2]+1),3)
  }
  
  if(tolower(method) %in% c("z-unpooled","unpooled","wald")){
    TX<-matrix(c(x,y,(p1-p2)/sqrt(p2*(1-p2)/Ns[2]+(p1)*(1-p1)/Ns[1])),(Ns[1]+1)*(Ns[2]+1),3)
  }
  
  if(tolower(method) %in% c("santner and snell","santner","snell")){
    TX<-matrix(c(x,y,(p1-p2)),(Ns[1]+1)*(Ns[2]+1),3)
  }
  
  #Observed Statistic:
  if(!(tolower(method) %in% c("csm","csm modified","csm approximate"))){
    TX[which(is.na(TX))]=0
    TXO<-TX[(Ns[2]+1)*data[1,1]+(data[1,2]+1),3]
  }
  
  if(tolower(method)=="csm approximate"){
    twosidedLess=(data[1,1]/Ns[1] <= data[1,2]/Ns[2])
    Tbls={}
    if (tolower(alternative)=="two.sided"){
      TX<-matrix(c(x,y,apply(matrix(c(x,Ns[1]-x,y,Ns[2]-y),(Ns[1]+1)*(Ns[2]+1),4),
                             1,FUN=function(tbls){maxPvalue(matrix(tbls,2,2))$pvalue})),(Ns[1]+1)*(Ns[2]+1),3)
    } else {
      otherTbls={}
      for(a in 0:Ns[1]){
        for(c in 0:Ns[2]){
          if(tolower(alternative)=="less" | (alternative=="two.sided" & twosidedLess)){
            if((a<data[1,1] & c>=data[1,2])|(a<=data[1,1] & c>data[1,2])){Tbls<-rbind(Tbls,c(a,c))
            } else if ((a<data[1,1] & c<data[1,2])|(a>data[1,1] & c>data[1,2])){otherTbls<-rbind(otherTbls,c(a,c))}
          } else {
            if((a>data[1,1] & c<=data[1,2])|(a>=data[1,1] & c<data[1,2])){Tbls<-rbind(Tbls,c(a,c))
            } else if ((a>data[1,1] & c>data[1,2])|(a<data[1,1] & c<data[1,2])){otherTbls<-rbind(otherTbls,c(a,c))}
          }
        }}
      otherTbls<-rbind(otherTbls,data[1,])
      TX<-matrix(c(otherTbls,apply(otherTbls,1,FUN=function(tbls){maxPvalue(matrix(tbls,1,2))$pvalue})),dim(otherTbls)[1],3)
    }
    
    if(!is.matrix(TX)){TX<-matrix(TX,ncol=3)}
    TXO<-max(TX[apply(matrix(TX[,1:2],ncol=2), 1, function(x) all(x == data[1,])),3])
    Tbls<-rbind(Tbls,TX[TX[,3] <= TXO,1:2])
  }
  
  #CSM test:
  if(tolower(method) %in% c("csm","csm modified")){
    
    #When the test is two-sided, do a one-sided test then add the other more
    #extreme tables at the end.  Have to choose which one-sided test to do though
    twosidedLess=(data[1,1]/Ns[1] <= data[1,2]/Ns[2])
    
    if(tolower(method)=="csm"){
      Tbls={}
      if(tolower(alternative)=="less" | (alternative=="two.sided" & twosidedLess)){Tbls<-matrix(c(0,Ns[2]),1,2)
      } else { Tbls<-matrix(c(Ns[1],0),1,2)}
      
      #AC represents the pair(s) of (a,c) being checked
      if(tolower(alternative)=="less" | (alternative=="two.sided" & twosidedLess)){
        AC<-matrix(c(0,Tbls[1,2]-1,Tbls[1,1]+1,Ns[2]),2,2,byrow=T)
        if(AC[2,1]>Ns[1]){AC=as.matrix(AC[-2,])}
        if(AC[1,2]<0){AC=as.matrix(AC[-1,])}
      } else {
        AC<-matrix(c(Ns[1],Tbls[1,2]+1,Tbls[1,1]-1,0),2,2,byrow=T)
        if(AC[2,1]<0){AC=as.matrix(AC[-2,])}
        if(AC[1,2]>Ns[2]){AC=as.matrix(AC[-1,])}
      }
    } else {
      Tbls={}
      for(a in 0:Ns[1]){
        for(c in 0:Ns[2]){
          if(tolower(alternative)=="less" | (alternative=="two.sided" & twosidedLess)){
            if((a<data[1,1] & c>=data[1,2])|(a<=data[1,1] & c>data[1,2])){Tbls<-rbind(Tbls,c(a,c))}
          } else {
            if((a>data[1,1] & c<=data[1,2])|(a>=data[1,1] & c<data[1,2])){Tbls<-rbind(Tbls,c(a,c))}
          }
        }}
      
      #AC represents the pair(s) of (a,c) being checked
      if(tolower(alternative)=="less" | (alternative=="two.sided" & twosidedLess)){
        AC<-matrix(c(0,data[1,2]-1,data[1,],data[1,1]+1,Ns[2]),3,2,byrow=T)
        if(AC[1,2]<0){AC=as.matrix(AC[-1,])}
        if(AC[dim(AC)[1],1]>Ns[1]){AC=as.matrix(AC[-dim(AC)[1],])}
      } else {
        AC<-matrix(c(Ns[1],data[1,2]+1,data[1,],data[1,1]-1,0),3,2,byrow=T)
        if(AC[1,2]>Ns[2]){AC=as.matrix(AC[-1,])}
        if(AC[dim(AC)[1],1]<0){AC=as.matrix(AC[-dim(AC)[1],])}
      }
    }
    TXO=NA
    
    if(is.null(Tbls)){Tbls<-matrix(data[1,],ncol=2)}
    
    if(tolower(alternative)=="two.sided"){
      Tbls<-rbind(Tbls,matrix(c(Ns[1]-Tbls[,1],Ns[2]-Tbls[,2]),dim(Tbls)[1],2))
      if(any(duplicated(Tbls))){Tbls<-matrix(Tbls[-which(duplicated(Tbls)),],ncol=2)}
    }
    
    while(sum(apply(Tbls, 1, function(x) all(x == data[1,])))==FALSE){
      
      #Calculate the possible more extreme test statistic:
      CcondAC<-rep(0,dim(AC)[1])
      for(j in 1:dim(AC)[1]){
        if(tolower(alternative)=='two.sided'){CcondAC[j]<-maxPvalue(rbind(Tbls,AC[j,],c(Ns[1]-AC[j,1],Ns[2]-AC[j,2])))$pvalue
        } else {CcondAC[j]<-maxPvalue(rbind(Tbls,AC[j,]))$pvalue}
      }
      
      addRow<-which(round(CcondAC,digits=10)==min(round(CcondAC,digits=10)))
      tempAdd<-matrix(AC[addRow,],length(addRow),2)
      
      #If table to be added is already in Tbls, ignore
      add<-NULL
      for(i in 1:dim(tempAdd)[1]){
        if(dim(merge(matrix(tempAdd[i,],ncol=2),Tbls))[1]==0){add<-rbind(add,tempAdd[i,])}
      }
      if(!is.null(add)){
        if(tolower(alternative)=="two.sided"){
          Tbls<-rbind(Tbls,add,matrix(c(Ns[1]-add[,1],Ns[2]-add[,2]),length(addRow),2))
          if(any(duplicated(Tbls))){Tbls<-matrix(Tbls[-which(duplicated(Tbls)),],ncol=2)}
        } else {
          Tbls<-rbind(Tbls,add)
        }
      }
      AC<-matrix(AC[-addRow,],ncol=2)
      
      if(alternative=="less" | (alternative=="two.sided" & twosidedLess)){
        newAdd<-matrix(c(add[,1],add[,1]+1,add[,2]-1,add[,2]),ncol=2)
        #It is possible that adding one success or one failure changes the sign, which shouldn't be considered 
        newAdd<-matrix(newAdd[apply(newAdd,1,function(x){(x[1]/Ns[1] <= x[2]/Ns[2])}),],ncol=2)
      } else {
        newAdd<-matrix(c(add[,1],add[,1]-1,add[,2]+1,add[,2]),ncol=2)
        newAdd<-matrix(newAdd[apply(newAdd,1,function(x){(x[1]/Ns[1] >= x[2]/Ns[2])}),],ncol=2)
      }
      AC<-rbind(AC,newAdd[!is.na(newAdd[,1]),])
      if(any(duplicated(AC))){AC<-matrix(AC[-which(duplicated(AC)),],ncol=2)}
      
      if(alternative=="less" | (alternative=="two.sided" & twosidedLess)){
        AC<-matrix(AC[order(AC[,1],-AC[,2]),],ncol=2)
        AC<-matrix(AC[!duplicated(AC[,1]),],ncol=2)
        AC<-matrix(AC[order(AC[,2],AC[,1]),],ncol=2)
        AC<-matrix(AC[!duplicated(AC[,2]),],ncol=2)
      } else {
        AC<-matrix(AC[order(AC[,1],AC[,2]),],ncol=2)
        AC<-matrix(AC[!duplicated(AC[,1]),],ncol=2)
        AC<-matrix(AC[order(AC[,2],-AC[,1]),],ncol=2)
        AC<-matrix(AC[!duplicated(AC[,2]),],ncol=2)
      }
      AC<-matrix(AC[AC[,1]>=0 & AC[,1]<=Ns[1] & AC[,2]>=0 & AC[,2]<=Ns[2],],ncol=2)
    }
    
  }
  
  #Find tables that have a test statistic as or more extreme than the observed statistic:
  if(!(tolower(method) %in% c("boschloo","csm","csm modified","csm approximate"))){
    if(tolower(alternative)=="greater"){Tbls<-matrix(TX[which(TX[,3]>=TXO),],ncol=3)
    } else if(tolower(alternative)=="less"){Tbls<-matrix(TX[which(TX[,3]<=TXO),],ncol=3)
    } else if(tolower(alternative)=="two.sided"){Tbls<-matrix(TX[which(abs(TX[,3])>=abs(TXO)),],ncol=3)}
  } else if (tolower(method) %in% c("boschloo")){Tbls<-matrix(TX[which(TX[,3]<=TXO),],ncol=3)}
  
  #Search for the maximum p-value:
  if(any(duplicated(Tbls))){stop("TABLES SHOULD NOT BE DUPLICATED: CHECK CODE")}
  maxP<-maxPvalue(Tbls)
  prob<-maxP$prob
  pvalue<-maxP$pvalue
  np<-maxP$np
  
  #Refine the p-value using the optimise function
  if(ref.pvalue==TRUE){
    refPvalue<-rep(0,length(np))
    refNp<-rep(0,length(np))
    for(i in 1:length(np)){
      ref<-optimise(f=function(p){sum(dbinom(Tbls[,1],Ns[1],p)*dbinom(Tbls[,2],Ns[2],p))},
                    interval=c(max(0.00001,np[i]-1/npNumbers),min(0.99999,np[i]+1/npNumbers)), maximum=TRUE)
      refPvalue[i]<-ref$objective+beta
      refNp[i]<-ref$maximum
    }
    np<-refNp[refPvalue==max(refPvalue)]
    pvalue<-max(refPvalue[refPvalue==max(refPvalue)])
  }
  
  #Note: if beta is large, can have a p-value greater than 1.  Cap at 1
  pvalue[pvalue>1]<-1
  
  #Plot p-value vs np
  if(to.plot==TRUE){
    plot(int,prob+beta,xlim=c(floor(min(int)*10)/10,ceiling(max(int)*10)/10),
         ylim=c(0,max(pvalue)),xlab="np",ylab="P-value",main="P-value as a function of the nuisance parameter")
    points(np,rep(pvalue,length(np)),col="red",pch=21,bg="red")
  }
  
  list<-list(method=method,p.value=max(pvalue),test.statistic=TXO,np=np,np.range=c(min(int),max(int)))
  return(list)
}
