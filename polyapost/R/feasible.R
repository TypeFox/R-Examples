#this function finds a feasible solution of A1x=b1,A2x<=b2,A3%*x>=b3,x>=eps or returns a negative vector when there is no feasible solution
#A1 is m1xq,A2 is m2xq,A3 is m3xq, b1 is m1x1 , b2 is m2x1,b3 is m3xq, eps is a scalar between 0 and 1 and has to be bigger than min(b3)

feasible<-function(A1,A2=NULL,A3=NULL,b1,b2=NULL,b3=NULL,eps)
{
  if (is.matrix(A1)){
      m1 <- nrow(A1)
      q<-ncol(A1)
  }
  else {
      m1 <- 1
      q<-length(A1)
  }
  if (!is.null(A2))
      if (is.matrix(A2))
          m2 <- nrow(A2)
      else m2 <- 1
  else m2 <- 0
  if (!is.null(A3))
      if (is.matrix(A3))
          m3 <- nrow(A3)
      else m3 <- 1
  else m3 <- 0

  if(m2+m3==0){
     A<-A1
     b<-b1
     out<-feasible1(A,b,rep(eps,q))
  }
  else {
     if (m2==0){
       A<-cbind(rbind(A1,A3),rbind(matrix(0,m1,m3),diag(-1,m3)))
       b<-c(b1,b3)
       out<-feasible1(A,b,rep(eps,q+m3))[1:q]
     }
     else {
        if(m3==0){
	  A<-cbind(rbind(A1,A2),rbind(matrix(0,m1,m2),diag(1,m2)))
	  b<-c(b1,b2)
	  out<-feasible1(A,b,rep(eps,q+m2))[1:q]

	}
	else{
	  A<-cbind(rbind(A1,A2,A3),rbind(matrix(0,m1,m2+m3),cbind(diag(1,m2),
	     matrix(0,m2,m3)),cbind(matrix(0,m3,m2),diag(-1,m3))))
          b<-c(b1,b2,b3)
        out<-feasible1(A,b,rep(eps,q+m2+m3))[1:q]
	}
     }
  }
  return(out)
}


