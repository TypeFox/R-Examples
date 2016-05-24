burtTable <- function(data) 
{
# burtTable makes a Burt matrix out of a data-frame 

nvar <- dim(data)[2]
ff <- rep(1, dim(data)[1])

ncat <- sapply(1:nvar,function(j) length(table(data[,j])))
tcat <- sum(ncat)
first <- cumsum(c(1,ncat))[1:nvar]
burt <- matrix(0,tcat,tcat)
vnames <- colnames(data)
namevec <- list()
for (k in 1:nvar) {
  namevec[[k]] <- paste(vnames[k], levels(data[,k]), sep=".")
}
rcnames <- unlist(namevec)
colnames(burt) <- rownames(burt) <- rcnames

for (i in 1:nvar) {
	ii <- first[i]:(first[i]+ncat[i]-1) 
  dd <- as.factor(data[,i])
	gi <- ifelse(outer(dd,levels(dd),"=="),1,0)
	if (i < nvar) {
		for (j in (i+1):nvar) {
			jj <- first[j]:(first[j]+ncat[j]-1)
      dd <- as.factor(data[,j])
			gj<-ifelse(outer(dd,levels(dd),"=="),1,0)
			burt[jj,ii]<-t(burt[ii,jj]<-crossprod(gi,gj*ff))
			}
		}
	burt[ii,ii]<-crossprod(gi,gi*ff)
	}

return(burt)
}