ctmat <-
function(Gmat,y,tr){
#creates a matrix (tmat) with final information of each terminal node of the tree (each column of Gmat)
#tmat= I*7 matrix
#Gmat = nodeindicator matrix
#y = outcome variable, column vector
#tr = treatmentvariable with two values (1: T=1; 2: T=2)
## cardinalities t1, cardinalities t2, and mean  and var y|t=1, mean and var y|t=2  
#each row of tmat gives this information for each column of Gmat
#thus number of rows of pmat corresponds to number of columns of Gmat
        rownum<-ncol(Gmat)
        t2mat<-matrix(0,ncol=2,nrow=rownum)
        t1mat<-matrix(0,ncol=2,nrow=rownum)
        tmat<-matrix(0,ncol=2,nrow=rownum)
        ##first column displays mean of y
        dat<-data.frame(cbind(y=y,tr=tr) )
      
        ##second column displays var of y
               t1mat[,1]<-sapply(1:rownum,function(kk,Gmat,y,tr){ifelse(sum(Gmat[,kk]==1&tr==1)==0,NA,mean(y[Gmat[,kk]==1&tr==1]))},Gmat=Gmat,y=y,tr=tr)
        t1mat[,2]<-sapply(1:rownum,function(kk,Gmat,y,tr){ifelse(sum(Gmat[,kk]==1&tr==1)==0,NA,sqrt(var(y[Gmat[,kk]==1&tr==1]))) },Gmat=Gmat,y=y,tr=tr)
         t2mat[,1]<-sapply(1:rownum,function(kk,Gmat,y,tr){ifelse(sum(Gmat[,kk]==1&tr==2)==0,NA,mean(y[Gmat[,kk]==1&tr==2])) },Gmat=Gmat,y=y,tr=tr)
        t2mat[,2]<-sapply(1:rownum,function(kk,Gmat,y,tr){ifelse(sum(Gmat[,kk]==1&tr==2)==0,NA,sqrt(var(y[Gmat[,kk]==1&tr==2]))) },Gmat=Gmat,y=y,tr=tr)
        tmat<-cbind(apply(Gmat[tr==1,],2,sum),t1mat,apply(Gmat[tr==2,],2,sum),t2mat)
        es<-sapply(1:rownum,function(kk,tmat){ifelse(is.na(sum(tmat[kk,c(2:3,5:6)])),NA,computeD(tmat[kk,1],tmat[kk,2],tmat[kk,3],tmat[kk,4],tmat[kk,5],tmat[kk,6])$dval) },tmat=tmat)
        se<-sapply(1:rownum,function(kk,tmat){ifelse(is.na(sum(tmat[kk,c(2:3,5:6)])),NA,computeD(tmat[kk,1],tmat[kk,2],tmat[kk,3],tmat[kk,4],tmat[kk,5],tmat[kk,6])$se) },tmat=tmat)
        #pval<-sapply(1:rownum,function(kk,tmat,Gmat,y,tr){ifelse(is.na(sum(tmat[kk,c(3:6)])),NA,t.test(y[Gmat[,kk]==1&tr==1],y[Gmat[,kk]==1&tr==0])$p.value) },tmat=tmat,Gmat=Gmat,y=y,tr=tr)
        tmat<-cbind(tmat,es,se  )
        colnames(tmat)<-c("nt1","meant1","sdt1","nt2","meant2","sdt2","d","se")
         return(as.matrix(tmat))}
