#this function transforms the constraints A1x=b1,A2x<=b2,A3%*x>=b3,x>=0 into #A1x=b1,A4x<=b4
constrmat<-function(A1,A2=NULL,A3=NULL,b2=NULL,b3=NULL)
{
  if (is.matrix(A1)){
      m1 <- nrow(A1)
      q<-ncol(A1)
  }
  else{
      m1 <- 1
      m1<-length(A1)
  }
  u<-rep(0,q)
  U<-diag(-1,q)
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
     A4<-U
     b4<-u
  }
  else {
     if (m2==0){
       A4<-rbind(-A3,U)
       b4<-c(-b3,u)
     }
     else {
        if(m3==0){
           A4<-rbind(A2,U)
	   b4<-c(b2,u)
        }
	else{
	  A4<-rbind(A2,-A3,U)
	  b4<-c(b2,-b3,u)
	}
     }
  }
  list(A1,A4,b4)
}
