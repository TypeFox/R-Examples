simulpln <- function(n,nitem,ncat,alphas,betas) {
  if(length(betas)!=nitem){stop("Length of betas does not match number of items")}
  if(length(alphas)!=nitem*(ncat-1)){stop("Length of alphas does not match nitem*(ncat-1)")}
    
  ## Size may differ from what C wants to return
  ## But, this is the max size that may be returned
  datvec<-rep(0,n*(nitem+1))
    
  tem <- .C("Rsimulpln", as.integer(nitem), 
     as.integer(ncat),as.integer(n),
     as.double(alphas),as.double(betas),
     datvec=as.double(datvec))
    
  ##nrec<-length(tem$datvec)/(nitems+1)
  ##datfr<-matrix(tem$datvec,nrec,nitem+1)
  datfr<-matrix(tem$datvec, ncol=nitem+1, byrow=TRUE)
    
  ## trim rows with 0 frequencies
  datfr<-datfr[which(datfr[,ncol(datfr)]>0),]    
    
  datfr
}