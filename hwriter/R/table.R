## private
## process .row or .col argument on its corresponding rcdata matrix
## location rule:
## - if rowcol named, expand it to dim(data)[kaxis], with first '' to headers
## - if rowcol unnamed of size 1, apply to header
## fill rule:
## - start at first data (not full) cell, if possible
## - recycle is TRUE by default
hprocessRowCol=function(rowcol,rcdata,data,axis,ddim,recycle=TRUE) {
  if (!is.null(rowcol)) {
    dim=dim(data)
    if (axis==2) ndata=rownames(data)
    else ndata=colnames(data)
    kaxis=3-axis
    
    ## named case
    nrowcol=names(rowcol)
    if (!is.null(nrowcol)) {
      if (is.null(ndata)) ndata=c('',rep(NA,dim[kaxis]-1))
      rowcol=rowcol[match(ndata,nrowcol)]
      if (is.list(rowcol)) rowcol=lapply(rowcol,function(z) if (is.null(z)) NA else z)
    } else {
      ## special rules to save room for the header:
      ## - if rowcol has more than one element
      ## - AND if one row/col has been added in the contraxis
      ## - AND if length rowcol is different from the contraxis
      if (length(rowcol)!=1 & dim[kaxis]!=ddim[kaxis] & length(rowcol)!=dim[kaxis]) rowcol=c(NA,rowcol)
    }
    
    for (i in 1:length(rowcol)) {
      rc=rowcol[[i]]
      if (!all(is.na(rc))) {
        if (recycle) rc=array(rc,dim=c(1,max(ddim[axis],length(rc))))
        drc=length(rc)
        j=ifelse(dim[axis]==drc,1,dim[axis]-ddim[axis]+1)
        if (axis==1) rcdata[j:(drc+j-1),i]=rc 
        else rcdata[i,j:(drc+j-1)]=rc
      }
    }
  }
  rcdata
}

## private
## reindex row or col named vector into a numeric fashion
hmapRowCol=function(x,rowcolnames) {
  if (is.null(x)) y=NULL
  else {
    if (is.null(names(x))) y=x
    else {
      if (is.null(rowcolnames)) y=NULL
      else {
        z=match(names(x),rowcolnames)
        y=rep(NA,length(rowcolnames))
        y[na.omit(z)]=x[which(!is.na(z))]
      }
    }
  }
  
  y
}

## handles row, col and table arguments
hwrite.table=function(data,page=NULL,...,table=TRUE,row.names=T,col.names=T,split.maxncol=NULL,split.maxnrow=NULL,col.width=NULL) { 
  ddim=dim(data)

  if (!table) {
    z=hwrite(as.vector(data),table=FALSE,...)
    dim(z)=ddim
    return(z)
  }
  
  ## process row and column names (add a col and a row, resp.)
  acol=row.names & !is.null(rownames(data))
  arow=col.names & !is.null(colnames(data))
  if (acol) data=cbind(rownames(data),data)
  if (arow) data=rbind(colnames(data),data)

  ## ddim is the data dim only (without extra col and/or row)
  ## dim  is the expanded data
  dim=dim(data)

  ## filters arguments
  args=list(...)
  nargs=names(args)
 
  ## table arguments
  itable=match(substr(nargs,1,6),'table.')
  iptable=match(nargs,c('cellspacing','cellpadding','width','border'))
  args.table=args[!is.na(itable)]
  names(args.table)=substr(names(args.table),7,nchar(names(args.table)))
  args.table=c(args.table,args[!is.na(iptable)])

  ## row arguments
  irow=match(substr(nargs,1,4),'row.')
  args.row=args[!is.na(irow)]
  names(args.row)=substr(names(args.row),5,nchar(names(args.row)))

  ## col arguments
  icol=match(substr(nargs,1,4),'col.')
  args.col=args[!is.na(icol)]
  names(args.col)=substr(names(args.col),5,nchar(names(args.col)))

  ## string arguments
  istring=match(nargs,c('name','heading','center','div','br'))
  args.string=args[!is.na(istring)]

  ## td arguments (remaining ones)
  args.td=args[is.na(itable) & is.na(irow) & is.na(icol) & is.na(istring) & is.na(iptable)]
  zargs=rep(list(NULL),length(args.row)+length(args.col))
  names(zargs)=c(names(args.row),names(args.col))
  args.td=c(args.td,zargs)
  
  ## expand if needed
  args.td=lapply(args.td,hexpand,dim,ddim)
  
  ## process .row and .cow arguments
  for (z in names(args.col)) args.td[[z]]=hprocessRowCol(args.col[[z]],args.td[[z]],data,1,ddim)
  for (z in names(args.row)) args.td[[z]]=hprocessRowCol(args.row[[z]],args.td[[z]],data,2,ddim)

  ## special case for col.width
  if (!is.null(col.width)) {
    width=array(NA,dim=dim)
    icol.width=hmapRowCol(col.width,colnames(data))
    width[1,1:length(icol.width)]= icol.width
    args.td=c(args.td,list(width=width))
  }

  ## process split.maxncol and split.maxnrow
  if (!is.null(split.maxncol) | !is.null(split.maxnrow)) {
    if (!is.null(split.maxncol) & !is.null(split.maxnrow)) stop('either \'split.maxnrow\' or \'split.maxncol\' must be NULL')
    
    ## split ancillary tables
    data=hsplitArray(data,maxnrow=split.maxnrow,maxncol=split.maxncol,preserve.size=T,output.list=F,arow=arow,acol=acol)
    for (z in names(args.td)) args.td[[z]]=hsplitArray(args.td[[z]],maxnrow=split.maxnrow,maxncol=split.maxncol,preserve.size=T,output.list=F,arow=arow,acol=acol)
  }
  
  do.call(hwriteRawTable,c(list(data,page=page,args.td=args.td,args.table=args.table,args.string=args.string)))
}

## private
## expands a to be of size ddb (if possible) and put in a matrix
## of size db, adding top/left NA rows/columns if needed
hexpand=function(a,db,ddb) {
  if (is.null(a)) a=NA
  if (is.null(dim(a))) a=array(a,dim=ddb)
  da=dim(a)
  b=array(NA,dim=db)
  i=ifelse(db[1]==da[1],1,db[1]-ddb[1]+1)
  j=ifelse(db[2]==da[2],1,db[2]-ddb[2]+1)
  b[i:(da[1]+i-1),j:(da[2]+j-1)]=a
  b
}

## private
hwriteRawTable=function(data,page=NULL,args.td=NULL,args.table=NULL,args.string=NULL) {
  ## default arguments
  if (is.null(args.table[['border']])) args.table[['border']]=1

  if (!is.matrix(data)) stop('\'data\' must be a matrix')
  dim=dim(data)
  data=as.vector(data)
  data[is.na(data)]='&nbsp;'
  
  ## process cell links
  data=hwrite(data,link=args.td$link,table=F)
  args.td$link=NULL
  
  ## process cells
  data=do.call(hmakeTag,c(list('td',data),args.td))
  
  ## process rows tr
  dim(data)=dim
  data=apply(data,1,function(z) paste(z,collapse=''))
  data=hmakeTag('tr',data,newline=T)

  ## process table
  data=paste(data,collapse='')
  str=do.call(hmakeTag,c(list('table',data,newline=T),args.table))
  
  ## final
  do.call(hwrite,c(list(str,page),args.string))
}

## private
## - split an array into a list of subarrays
## - rownames and colnames are ignored
## - arow and acol indicate that the first row (resp col) should be treated as a rowname (resp colname)
hsplitArray=function(data,maxnrow=0,maxncol=0,preserve.size=T,output.list=T,arow=F,acol=F) {
  if (is.null(data)) return(NULL)
  if (!is.matrix(data)) stop('\'data\' must be a matrix')
  if (is.null(maxnrow)) maxnrow=0
  if (is.null(maxncol)) maxncol=0
  if (maxnrow==0 && maxncol==0) stop('splitting must be done in one direction: either \'maxnrow\' or \'maxncol\' must be non null')
  if (maxnrow*maxncol!=0) stop('splitting cannot be done on both directions: either \'maxnrow\' or \'maxncol\' must be null')
  if (!output.list && !preserve.size) stop('outputting matrix is possible only if \'preserve.size\' is TRUE')

  ## preserve first row (resp col)
  if (arow) {
    zrow=data[1,]
    data=data[-1,]
  }
  if (acol) {
    zcol=data[,1]
    data=data[,-1]
  }
  
  nr=nrow(data)
  nc=ncol(data)
  if (output.list) out=list()
  else out=NULL
  
  ## maxncol splitting (horizontal splitting)
  if (maxncol>0) { 
    nc2=maxncol  
    z1=1
    repeat {
      if (z1>nc) break
      z2=z1+nc2-1
      if (z2>nc) z2=nc
      if (preserve.size) {
        napz=nc2-z2+z1-1
        zdata=matrix(c(data[,z1:z2],rep(NA,nr*napz)),nrow=nr,ncol=nc2)
        nap=rep(NA,napz)
      } else {
        zdata=matrix(data[,z1:z2],nrow=nr,ncol=1+z2-z1)
        nap=NULL
      }
      
      if (acol&arow) {
        zdata=cbind(zcol,zdata)
        zdata=rbind(c(zrow[c(1,1+(z1:z2))],nap),zdata)
      }
      else if (acol) zdata=cbind(zcol,zdata)
      else if (arow) zdata=rbind(c(zrow[z1:z2],nap),zdata)
     
      if (output.list) out=c(out,list(zdata))
      else out=rbind(out,zdata)
      z1=z2+1
    }
  }
  
  ## maxnrow splitting (vertical splitting)
  if (maxnrow>0) { 
    nr2=maxnrow
    z1=1
    repeat {
      if (z1>nr) break
      z2=z1+nr2-1
      if (z2>nr) z2=nr
      if (preserve.size) {
        napz=nr2-z2+z1-1
        zdata=matrix(c(t(data[z1:z2,]),rep(NA,nc*napz)),nrow=nr2,ncol=nc,byrow=T)
        nap=rep(NA,napz)
      } else {
        zdata=matrix(data[z1:z2,],nrow=1+z2-z1,ncol=nc)
        nap=NULL
      }
      if (acol) zdata=cbind(c(zcol[z1:z2],nap),zdata)
      if (arow) zdata=rbind(zrow,zdata)
      
      if (output.list) out=c(out,list(zdata))
      else out=cbind(out,zdata)
      z1=z2+1
    }
  }
  
  out
}
