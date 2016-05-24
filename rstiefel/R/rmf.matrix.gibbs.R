rmf.matrix.gibbs <-
function(M,X,rscol=NULL)
{
  #simulate from the matrix mf distribution as described in Hoff(2009) 
  #this is one Gibbs step, and must be used iteratively

  #The number of columns to replace should be small enough so that 
  #rmf.matrix will be quick, but not so small that you are not taking 
  #advantage of rmf.matrix. In particular, for square matrices you must 
  #be replacing two or more columns at a time, or else the chain is 
  #not irreducible. 

  if(is.null(rscol)){rscol<-max(2,min( round(log(dim(M)[1])),dim(M)[2])) }

  sM<-svd(M)
  H<-sM$u%*%diag(sM$d)
  Y<-X%*%sM$v

  m<-dim(H)[1] ; R<-dim(H)[2]
  for(iter in 1:round(R/rscol))
  {
    r<-sample(seq(1,R,length=R),rscol)
    N<-NullC(Y[,-r])
    y<-rmf.matrix(t(N)%*%H[,r])
    Y[,r]<- N%*%y
  }
  Y%*%t(sM$v)
}
