cdcor <-
function(x,y,z,width,index=1) {

  x<-as.matrix(x);y<-as.matrix(y);z<-as.matrix(z)
  dim_x<-dim(x); n<-dim_x[1];M<-dim_x[2];p<-1; q<-dim(y)[2];d<-dim(z)[2]
  k<-numeric(n*n)
  CDCOV<-numeric(n)
  iraCDCOV<-numeric(M)
  re <- .C("IterationcdCov", as.double((x)), as.double(t(y)), as.double(t(z)), as.integer(n), 
        as.integer(p), as.integer(q), as.integer(d), as.double(index),as.double(width), k=as.double(k),
        cd=as.double(CDCOV),as.integer(M),iracd=as.double(iraCDCOV)) 
  cdc <- list(cdcor=re$cd, mcdcor = re$iracd, width=width, index=index)
  return(cdc)
}
