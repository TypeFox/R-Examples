abc.a.coefs <-
function(indices, splitting=FALSE) {

  # im Unterschied zu a.coefs hier nur 1en in Matrix!!
  a.coefs.nominal <- function (p = NULL, ...)
      { if (p > 1)
        {
        a.coefs.mat <- matrix(nrow=p,ncol=2^p)  # sum(choose(p,1:p)) == 2^p
        for (i in 1:p) {
             a.coefs.mat[i,] <- rep(c(1,0), times = c(2^(p-i),2^(p-i)))
             }
        a.coefs.mat <- a.coefs.mat[,-c(2^p)]
        } else {
        a.coefs.mat <- matrix(1,ncol=1,nrow=1)        
        }
        return (a.coefs.mat) 
      }

  a.coefs.ordinal <- function (p = NULL, ...) 
      { if (p > 1)
        {

        if (p > 2)
          { a.coefs.mat <- diag(p)
            for (i in 1:(p-2))
               { m <- matrix(0,ncol=p, nrow=p-i)
                 for (j in 1:(i+1)) { m <- m + cbind (matrix(0,ncol=j-1,nrow=p-i), diag (p-i), matrix(0,ncol=i-j+1,nrow=p-i)) }
                 a.coefs.mat <- cbind(a.coefs.mat, t(m))
               }
            a.coefs.mat <- cbind(a.coefs.mat,1)
          }
        else { a.coefs.mat <- cbind (diag (2), c(1,1)) }
        } else {
        a.coefs.mat <- matrix(1,ncol=1,nrow=1)
        }

        return (a.coefs.mat) }
       
  index1 <- indices[[1]]
  index2 <- indices[[2]]

#  require(Matrix)
  A <- matrix(0,ncol=0,nrow=0)
  
  if (splitting==TRUE){  # index1, index2 dann nur 1 Eintrag, wie Bloecke von Laenge index1 zerlegt werden koennen
  for (i in 1:length(index1)) {
  if ( index2[i]<0 ) {  B <- a.coefs.nominal(p=index1[i])[,-1] 
                        B <- B[,which(B[1,]==1)]
                        A <- bdiag(A,B)} 
  if ( index2[i]==0) {  A <- bdiag(A,matrix(0,ncol=0,nrow=index1[i])) }
  if ( index2[i]>0 ) {  B <- a.coefs.ordinal(p=index1[i])[,-(sum(1:index1[i]))] 
                        B <- B[,which(B[1,]==1)]
                        A <- bdiag(A,B)} 
  } }  else {
  
  # splitting==FALSE
  for (i in 1:length(index1)) {  
  if ( index2[i]<0 ) {  A <- bdiag(A,a.coefs.nominal(p=index1[i])) }
  if ( index2[i]==0) {  A <- bdiag(A,matrix(0,ncol=0,nrow=index1[i])) }
  if ( index2[i]>0 ) {  A <- bdiag(A,a.coefs.ordinal(p=index1[i])) }
  } }

  a.coefs.mat<-as.matrix(A)
  return(a.coefs.mat)
}

