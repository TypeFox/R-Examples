`Penalty.matrix` <-function(m,order=2){
    ######################################################
    ## internal function that computes the differences
    Diff.matrix<-function(m,order=2){
        # internal function that computes the first order differences
        d.matrix<-function(m){

            A<-cbind(diag(m-1),rep(0,m-1))

            B<-cbind(rep(0,m-1),-1*diag(m-1))

            d<-A+B

        return(d)

        }

    D<-d.matrix(m)

  if (order>1){

    for (k in 2:order){

      D<-d.matrix(m-k+1)%*%D

    }

  }

  return(D)

}

########################################################################
 # m is a vector that contains the size of the blocks for each penalty term
  p<-length(m)
  
  start.block=cumsum(m) - m +1
  end.block=cumsum(m)
  P<-matrix(0,sum(m),sum(m))
  
  for (i in 1:p){
    D<-Diff.matrix(m[i],order=order)
    K<-t(D)%*%D
    P[start.block[i]:end.block[i],start.block[i]:end.block[i]]=K
}
  return(P)

}
