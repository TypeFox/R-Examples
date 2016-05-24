hwrite=function(x,page=NULL,...)
  UseMethod('hwrite')

hwrite.character=function(x,...)
  hwrite.vector(x,...)

hwrite.numeric=function(x,...)
  hwrite.vector(x,...)

hwrite.array=function(x,...)
  hwrite.table(x,...)

hwrite.matrix=function(x,...)
  hwrite.table(x,...)

hwrite.data.frame=function(x,...)
  hwrite.table(as.matrix(x),...)

## switch between hwriteString and hwrite.matrix
## redimension 'dim' and 'byrow' matrix orientation
hwrite.vector=function(data,page=NULL,...,table=NULL,names=NULL,byrow=NULL,dim=NULL) {
  ## default arguments
  if (is.null(table)) {
    if (length(data)<=1) table=FALSE
    else table=TRUE
  }
  if (is.null(names)) names=TRUE
  if (is.null(byrow)) byrow=FALSE
    
  if (table) {
    if (is.null(dim)) dim=c(1,length(data))
    datanames=names(data)
    data=matrix(data,nrow=dim[1],ncol=dim[2],byrow=byrow)
    mode(data)='character'
    ## preserve names, if possible
    if (names) {
      if (dim[1]==length(datanames)) rownames(data)=datanames
      if (dim[2]==length(datanames)) colnames(data)=datanames
    }
    hwrite.matrix(data,page=page,...)
  } else hwriteString(data,page=page,...)
}

## private
## final string writing function
hwriteString=function(txt,page=NULL,...,link=NULL,name=NULL,heading=NULL,center=NULL,br=NULL,div=NULL) {
  ## default arguments
  if (is.null(br)) br=FALSE
  if (is.null(center)) center=FALSE
  if (is.null(div)) div=FALSE
  args=list(...)

  ## box text with:
  ## - 'a'    if link is non-null or name is non-null
  ## - 'h*'   if heading is non-null
  ## - 'div'  if div is TRUE
  ## - 'span' if args are present
  ## - no box otherwise
  ##
  ## also: removes <a> tags if corresponding href and argument values are NA
  boxtag=NULL
  if (!is.null(link)) {
    args=c(args, list(href=link))
    boxtag = rep('a', length(link))
    boxtag[is.na(link)] = NA
  }
  else if (!is.null(name)) {
    args=c(args, list(name=name))
    boxtag = rep('a', length(name))
    boxtag[is.na(name)] = NA
  }
  else if (!is.null(heading)) boxtag=paste('h',heading,sep='')
  else if (div) boxtag='div'
  else if (length(args)>0) boxtag='span'

  ## box text
  if (!is.null(boxtag)) txt=do.call(hmakeTag,c(list(boxtag,txt),args))
  
  ## center
  if (center) txt=hmakeTag('center',txt)
  
  ## line break
  if (br) txt=paste(txt,'<br/>\n',sep='')

   ## final output
  if (is.null(page)) txt
  else if (is.character(page)) {
    p=openPage(page)
    cat(txt,file=p)
    closePage(p)
    invisible(txt)
  } else invisible(cat(txt,file=page)) 
}

hwriteImage=function(image.url,page=NULL,...,image.border=0,width=NULL,height=NULL,capture=FALSE) {
  ## take a snapshot of the current device ?
  if (capture) {
    if (is.null(width)) width=400
    if (is.null(height)) height=400
    dev.print(png,width=width,height=height,image.url)
  }
  str=hmakeTag('img',border=image.border,src=image.url,alt=image.url,width=width,height=height)
  
  ## final
  hwrite(str,page,...)
}

resync=function() {
  try(detach('package:hwriter'),silent=TRUE)
  hwrite=NULL
  source('R/hwriter.R')
  source('R/page.R')
  source('R/table.R')
  source('R/example.R')
  library(hwriter)
}

hmakeTag = function(tag, data=NULL, ..., newline=FALSE) {
  attrs = list(...)

  ## dim is the output dim of the result
  dim = dim(tag)
  if (!is.null(dim(data))) dim = dim(data)

  if (is.null(data)) data = ''
  na = length(attrs)

  ## attributes grid
  xattrs = NULL
  if (na>0) {
    namax = max(sapply(attrs, length))
    n = max(c(length(tag), length(data), namax))
    xattrs = matrix('', nrow=n, ncol=na)
    nattrs = names(attrs)
    for (i in 1:na) {
      z = attrs[[i]]
      if (!is.null(z)) {
        fna = !is.na(z)
        xattrs[fna,i] = paste(' ',nattrs[i], '=\"', z[fna], '\"', sep='')
        if (!is.null(dim(z))) dim = dim(z)
      }
    }
    xattrs = apply(xattrs, 1, paste, collapse='')
  }
  
  if (newline) nl = '\n' else nl = NULL
  res = paste('<', tag, xattrs, '>', nl, data, '</', tag, '>', nl, sep='')
  natag = rep(is.na(tag), length(res)/length(tag))
  res[natag] = paste(rep('', length(tag)), rep('', length(xattrs)), data, sep='')[natag]
  if (!is.null(dim)) res = array(res, dim=dim)
  res
}
