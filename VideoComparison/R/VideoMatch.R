VideoMatch<-function(lh1,lh2,sc=1){
	aux1<-function(x,j){
	  return(x[j])
	}
#
	VideoVectorSearch<-function(x,l2){
  		rs<-lapply(l2,VideoDistance,x)
  		ev<-matrix(0,nrow=4,ncol=length(l2))
  		ul<-unlist(lapply(rs,aux1,1))
  		ev[1,]<-ul/max(ul)
  		ul<-unlist(lapply(rs,aux1,2))
  		ev[2,]<-ul/max(ul)
  		ul<-unlist(lapply(rs,aux1,3))
  		ev[3,]<-ul/max(ul)
  		ul<-unlist(lapply(rs,aux1,4))
  		ev[4,]<-1 - ul/max(ul)
  		mev<-apply(ev,2,sum)/4           
  		j<-which(mev==min(mev))  
  		return(list(idx=j,err=mev[j]))
	}
#
  if ( ! is.list(lh1) | ! is.list(lh2)) return(NULL);
  if ( length(lh1) > length(lh2)) {
    l1<-lh1
    l2<-lh2
    pos<-0
  } else {
    l1<-lh2
    l2<-lh1
    pos<-1    
  }
  rs<-lapply(l1,VideoVectorSearch,l2)
  unord<-unlist(lapply(rs,function(x){return(x$idx)}))
  unerr<-unlist(lapply(rs,function(x){return(x$err)}))
  names(unerr)<-1:length(unerr)
  vl<-tapply(unerr,unord,function(x){return(names(which.min(x)))})
  unord<-sum(diff(as.numeric(vl))<0)
  err<-mean(unerr[as.numeric(vl)])
  fct<-sc
  if (unord > 0 ) fct<-fct*0.5
  return((1-err)*fct)
}

