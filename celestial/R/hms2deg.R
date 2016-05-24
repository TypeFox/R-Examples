hms2deg <-
function(h=0,m=0,s=0,sep=':'){
  if(length(dim(h))==2){
    if(dim(h)[2]==3){
      if(is.character(h[1,1]) & missing(m) & missing(s)){
        m=as.numeric(h[,2])
        s=as.numeric(h[,3])
        h=h[,1]
        signlogic=grepl('-',h)
        sign=rep(1,length(h))
        sign[signlogic]=-1
        h=abs(as.numeric(h))
      }else{
        m=h[,2]
        s=h[,3]
        h=h[,1]
      }
    }else{stop("d has wrong dimension, should be a Nx3 table/matrix")}
  }
  if(is.character(h[1]) & missing(m) & missing(s)){
    if(sep!='HMS' & sep!='hms'){
      temp=strsplit(h,split=sep,fixed=TRUE)
      split=as.numeric(unlist(temp))
      skip=length(temp[[1]])
    }
    if(sep=='HMS'){split=unlist(strsplit(h,split='H',fixed=TRUE));skip=2}
    if(sep=='hms'){split=unlist(strsplit(h,split='h',fixed=TRUE));skip=2}
    nsplit=length(split)/skip
    h=as.numeric(split[seq(1,(nsplit-1)*skip+1,by=skip)])
    
    if(sep!='HMS' & sep!='hms'){
      if(skip>=2){m=as.numeric(split[seq(2,(nsplit-1)*skip+2,by=skip)])}
      if(skip>=3){s=as.numeric(split[seq(3,(nsplit-1)*skip+3,by=skip)])}
    }
    if(sep=='HMS'){
    split=unlist(strsplit(split[seq(2,(nsplit-1)*skip+2,by=skip)],split='M',fixed=TRUE))
    m=as.numeric(split[seq(1,(nsplit-1)*skip+1,by=skip)])
    s=as.numeric(unlist(strsplit(split[seq(2,(nsplit-1)*skip+2,by=skip)],split='S',fixed=TRUE)))
    }
    if(sep=='hms'){
    split=unlist(strsplit(split[seq(2,(nsplit-1)*skip+2,by=skip)],split='m',fixed=TRUE))
    m=as.numeric(split[seq(1,(nsplit-1)*skip+1,by=skip)])
    s=as.numeric(unlist(strsplit(split[seq(2,(nsplit-1)*skip+2,by=skip)],split='s',fixed=TRUE)))
    }
  }
if(any(h<0 | h>24)){stop('All h values should be 0<=h<=24')}
if(any(m<0 | m>=60)){stop('All m values should be 0<=m<60')}
if(any(s<0 | s>=60)){stop('All s values should be 0<=s<60')}
    if(is.matrix(h) || is.data.frame(h)){
		if(ncol(h) == 1){h=h[,1]}
        else	if(ncol(h) == 2){m = h[, 2];h = h[, 1]}
        else	if(ncol(h) == 3){s = h[, 3];m = h[, 2];h = h[, 1]}
        }
    H=floor(as.numeric(h))
    M=floor(as.numeric(m))
    S=as.numeric(s)
    totaldeg=(H*15) + (M*15/60) + (S*15/3600)
    return(totaldeg)
}
