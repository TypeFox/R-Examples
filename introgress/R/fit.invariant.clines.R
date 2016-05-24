fit.invariant.clines <-
function(sam=NULL,n.ind=NULL,locustype=NULL){
  if (locustype=="C" | locustype=="c"){
    if(sum(sam==2,na.rm=TRUE)>=1){
      AA.fitted<-rep(1,n.ind)
      Aa.fitted<-rep(0,n.ind)
      aa.fitted<-rep(0,n.ind)
    }
    else if(sum(sam==1,na.rm=TRUE)>=1){
      AA.fitted<-rep(0,n.ind)
      Aa.fitted<-rep(1,n.ind)
      aa.fitted<-rep(0,n.ind)
    }
    else if(sum(sam==0,na.rm=TRUE)>=1){
      AA.fitted<-rep(0,n.ind)
      Aa.fitted<-rep(0,n.ind)
      aa.fitted<-rep(1,n.ind)
    }
    fitted.mat<-rbind(AA.fitted,Aa.fitted,aa.fitted)
    return(fitted.mat)
  }
  else {
    if(sum(sam==1,na.rm=TRUE)>=1){
      AA.fitted<-rep(1,n.ind)
      aa.fitted<-rep(0,n.ind)
    }
    else if(sum(sam==0,na.rm=TRUE)>=1){
      AA.fitted<-rep(0,n.ind)
      aa.fitted<-rep(1,n.ind)
    }
    fitted.mat<-rbind(AA.fitted,aa.fitted)
    return(fitted.mat)
  }
}

