mal.eq <- function(S,P,N){
  phi<-diag(0/N)
  Pt<-t(P)
  x<-0
  repeat{
    x<-x+1
    S1<-mtx.exp(S,x)
    P1<-mtx.exp(P,x)
    Pt1<-mtx.exp(Pt,x)
    D<-(1-phi)/(2*N)
    D<-diag(D)
    D<-diag(D) ## everything till here is similar to a normal phi calculation
    toll<-phi ## I use toll as a comparison mark. toll is phi a n-1 cycles
    toll1<-signif(toll,6) ## optional. I set the number of significant digits to 6
    phi<-phi+(S1%*%Pt1%*%D%*%P1%*%S1) ## that's phi at n cycles
    phi1<-signif(phi,6) ## optional. As for toll
    if (identical(toll1,phi1)){ ## logical condition. If toll (that is, phi for n-1) and phi are identical
      return(x-1) ## return the value of n-1
      break ## and stop, because the Malecot model has reached its asymptot
    }
  }
}
