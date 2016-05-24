dms2deg <-
function(d=0,m=0,s=0,sign='d',sep=':'){
  if(length(dim(d))==2){
    if(dim(d)[2]==3){
      if(is.character(d[1,1]) & missing(m) & missing(s)){
        m=as.numeric(d[,2])
        s=as.numeric(d[,3])
        d=d[,1]
        signlogic=grepl('-',d)
        sign=rep(1,length(d))
        sign[signlogic]=-1
        d=abs(as.numeric(d))
      }else{
        m=d[,2]
        s=d[,3]
        d=d[,1]
      }
    }else{stop("d has wrong dimension, should be a Nx3 table/matrix")}
  }
  if(is.character(d[1]) & missing(m) & missing(s)){
    if(sep!='DMS' & sep!='dms'){
      temp=strsplit(d,split=sep,fixed=TRUE)
      split=as.numeric(unlist(temp))
      skip=length(temp[[1]])
    }
    if(sep=='DMS'){split=unlist(strsplit(d,split='D',fixed=TRUE));skip=2}
    if(sep=='dms'){split=unlist(strsplit(d,split='d',fixed=TRUE));skip=2}
    nsplit=length(split)/skip
    d=split[seq(1,(nsplit-1)*skip+1,by=skip)]
    signlogic=grepl('-',d)
    sign=rep(1,nsplit)
    sign[signlogic]=-1
    d=abs(as.numeric(d))
    if(sep!='DMS' & sep!='dms'){
      if(skip>=2){m=as.numeric(split[seq(2,(nsplit-1)*skip+2,by=skip)])}
      if(skip>=3){s=as.numeric(split[seq(3,(nsplit-1)*skip+3,by=skip)])}
    }
    if(sep=='DMS'){
    split=unlist(strsplit(split[seq(2,(nsplit-1)*skip+2,by=skip)],split='M',fixed=TRUE))
    m=as.numeric(split[seq(1,(nsplit-1)*skip+1,by=skip)])
    s=as.numeric(unlist(strsplit(split[seq(2,(nsplit-1)*skip+2,by=skip)],split='S',fixed=TRUE)))
    }
    if(sep=='dms'){
    split=unlist(strsplit(split[seq(2,(nsplit-1)*skip+2,by=skip)],split='m',fixed=TRUE))
    m=as.numeric(split[seq(1,(nsplit-1)*skip+1,by=skip)])
    s=as.numeric(unlist(strsplit(split[seq(2,(nsplit-1)*skip+2,by=skip)],split='s',fixed=TRUE)))
    }
  }
if(any(d< -90 | d>90) & sign[1]=='d'){stop('All d values should be -90<=d<=90')}
if(any(d< 0 | d>90) & sign[1]!='d'){stop('Since sign is specified, all d values should be 0<=d<=90')}
if(any(m<0 | m>=60)){stop('All m values should be 0<=m<60')}
if(any(s<0 | s>=60)){stop('All s values should be 0<=s<60')}

if(sign[1]=='d' & any(d==0)){stop('Some d values are 0, and these have ambiguous signs, specify sign explicitly')}
if(sign[1]=='d'){sign=sign(d);d=abs(d)}
if(sign[1]!='d' & length(sign) != length(d)){stop('d and sign lengths do not match')}

    if(is.matrix(d) || is.data.frame(d)){
		if(ncol(d) == 1){d=d[,1]}
        else	if(ncol(d) == 2){m = d[, 2];d = d[, 1]}
        else	if(ncol(d) == 3){s = d[, 3];m = d[, 2];d = d[, 1]}
    }
    D=floor(as.numeric(d))
    M=floor(as.numeric(m))
    S=as.numeric(s)
    totaldeg=(D)*sign + (M/60)*sign + (S/3600)*sign
    return(totaldeg)
}
