fit.c.clines <-
function(reg.out=NULL, hi.index=NULL, sam=NULL, n.ind=NULL){
  if (length(coef(reg.out))>2){
    AA.slope<-coef(reg.out)[2,2]
    AA.int<-coef(reg.out)[2,1]
    Aa.slope<-coef(reg.out)[1,2]
    Aa.int<-coef(reg.out)[1,1]
    Hx<-exp(AA.slope*hi.index+AA.int) + exp(Aa.slope*hi.index+Aa.int)
    AA.fitted<-exp(AA.slope*hi.index+AA.int)/(1+Hx)
    Aa.fitted<-exp(Aa.slope*hi.index+Aa.int)/(1+Hx)
    aa.fitted<-1-(Aa.fitted+AA.fitted)
  }	
  else if (length(coef(reg.out))==2){
    if(sum(sam==2,na.rm=TRUE)>=1 & sum(sam==1,na.rm=TRUE)>=1){
      AA.slope<-coef(reg.out)[2]
      AA.int<-coef(reg.out)[1]
      Aa.slope<-NA
      Aa.int<-NA
      Hx<-exp(AA.slope*hi.index+AA.int)
      AA.fitted<-exp(AA.slope*hi.index+AA.int)/(1+Hx)
      Aa.fitted<-1-exp(AA.slope*hi.index+AA.int)/(1+Hx)
      aa.fitted<-rep(0,n.ind)
    }
    else if(sum(sam==2,na.rm=TRUE)>=1 & sum(sam==0,na.rm=TRUE)>=1){
      AA.slope<-coef(reg.out)[2]
      AA.int<-coef(reg.out)[1]
      Aa.slope<-NA
      Aa.int<-NA
      Hx<-exp(AA.slope*hi.index+AA.int)
      AA.fitted<-exp(AA.slope*hi.index+AA.int)/(1+Hx)
      Aa.fitted<-rep(0,n.ind)
      aa.fitted<-1-exp(AA.slope*hi.index+AA.int)/(1+Hx)
    }
    else if(sum(sam==1,na.rm=TRUE)>=1 & sum(sam==0,na.rm=TRUE)>=1){
      AA.slope<-NA
      AA.int<-NA
      Aa.slope<-coef(reg.out)[2]
      Aa.int<-coef(reg.out)[1]
      Hx<-exp(Aa.slope*hi.index+Aa.int)
      AA.fitted<-rep(0,n.ind)
      Aa.fitted<-exp(Aa.slope*hi.index+Aa.int)/(1+Hx)
      aa.fitted<-1-exp(Aa.slope*hi.index+Aa.int)/(1+Hx)
    }
  }
  fitted.mat<-rbind(AA.fitted,Aa.fitted,aa.fitted)
  return(fitted.mat)					
}

