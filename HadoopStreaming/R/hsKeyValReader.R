`hsKeyValReader` <-
function(file="",chunkSize=-1,skip=0, sep='\t',FUN=function(k,v) cat(paste(k,v,sep=': '),sep='\n')) {
  ## can val have separator characters?  
  if (skip>0) {
    junk = scan(file,what='character',quiet=TRUE,sep=sep,flush=TRUE,nlines=skip)
  }
  cols = list(keys='',vals='')
  repeat {
    a = scan(file,what=cols,quiet=TRUE,sep=sep,flush=TRUE,nlines=chunkSize)
    if (length(a$keys) ==0) break
    FUN(a$keys,a$vals)
  }
}

