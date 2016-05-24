cpmat <-
function(Gmat,y,tr,crit="es"){
#creates pmat = candidate parent nodes information matrix
#pmat= I * 3 matrix ; I = number of candidate parent nodes
#columns of pmat: per node: cardinality t1, cardinality t2, cohen's d or difference in absolute means
#Gmat = N*I matrix=indicator matrix of all candidate parent nodes
#y = outcome variable , column vector
#tr = treatment-variable with values 1 and 2, column vector
#crit="es": effect size criterion;  crit="dm": difference in means
        rownum<-ncol(Gmat)
        pmat<-matrix(0,ncol=3,nrow=rownum)
        nt2<-apply(Gmat[tr==2,],2,sum)
        nt1<-apply(Gmat[tr==1,],2,sum)
        meant2<-numeric(rownum)
        meant1<-numeric(rownum)
        #compute mean value of y for T=0 
        meant2<-sapply(1:rownum,function(kk,Gmat,y,tr){mean(y[Gmat[,kk]==1&tr==2]) },Gmat=Gmat,y=y,tr=tr)
        #print(meant0)
        #compute mean value of y for T=1
        meant1<-sapply(1:rownum,function(kk,Gmat,y,tr){mean(y[Gmat[,kk]==1&tr==1]) },Gmat=Gmat,y=y,tr=tr)
       # print(meant1)
        pmat<-as.matrix(cbind(nt1,nt2,meant1-meant2)) 
        if(crit=="es"){
         sigmap<-numeric(rownum)
        sigmap<-sapply(1:rownum,function(kk,Gmat,y,tr,nt1,nt2){
        ifelse(nt1[kk]<=1|nt2[kk]<=1,0,csigmap(nt1[kk],sqrt(var(y[Gmat[,kk]==1&tr==1])) ,
           nt2[kk],sqrt(var(y[Gmat[,kk]==1&tr==2]))) )},Gmat=Gmat,y=y,tr=tr,nt1=nt1,nt2=nt2)
        pmat[,3]<-pmat[,3]/sigmap 
             pmat[,3][pmat[,3]==Inf]<-0
             pmat[,3][pmat[,3]==-Inf]<-0}
          pmat[,3][is.na(pmat[,3])]<-0
          dimnames(pmat)<-NULL
        return(pmat)}
