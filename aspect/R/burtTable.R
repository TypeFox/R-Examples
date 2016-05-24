`burtTable` <-
function(data)
{
# creates burt matrix (each variable is assumed to be categorical)
  
  nvar <- dim(data)[2]
  ncat <- sapply(1:nvar,function(j) length(table(data[,j])))
  tcat <- sum(ncat)
  first<-cumsum(c(1,ncat))[1:nvar]
  burt <- matrix(0,tcat,tcat)
  for (i in 1:nvar) {
    ii<-first[i]:(first[i]+ncat[i]-1)
    if (i < nvar) {
      for (j in (i+1):nvar) {
	jj<-first[j]:(first[j]+ncat[j]-1)
	burt[jj,ii]<-t(burt[ii,jj]<-as.matrix(table(data[,i],data[,j])))
      }
    }
    burt[ii,ii]<-as.matrix(table(data[,i],data[,i]))
  }
return(burt)
}

