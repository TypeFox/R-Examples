cpos <- function(str,sub,start=1)
{### find the first position of string sub in string str, starting from position start
  lstr  <- nchar(str)
  lsub1 <- nchar(sub)-1
  if (start+lsub1 > lstr) return(NA)
  else {
    str <- substring(str,start,lstr)
    str <- substring(str, 1:(lstr-lsub1), (1+lsub1):lstr)
    p <- ((start:lstr)[str==sub])[1]
    if (is.na(p>0)) return(NA)
    else return(p)
  }
}

cposR <- function(str, sub, restrict) {
  if(length(str)>1) stop('only works with a single string str')
  l.str   <- nchar(str)
  l.sub <- nchar(sub)
  if(l.sub > l.str) return(list(first=0,last=0))
  if(l.sub==l.str)  return(if(str==sub)list(first=1,last=l.str) else 
    list(first=0,last=0))

  is <- 1:(l.str-l.sub+1)
  ss <- substring(str, is, is+l.sub-1)
  k <- ss==sub
  if(!any(k)) return(list(first=NA,last=NA))
  k <- is[k]
  if(!missing(restrict)) k <- k[k>=restrict[1] & k<=restrict[2]]
  if(length(k)==0) return(list(first=NA,last=NA))
  list(first=k, last=k+l.sub-1)
}

