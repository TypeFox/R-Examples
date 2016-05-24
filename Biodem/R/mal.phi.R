mal.phi <- function(S,P,N,n){
  if (n < 1){
    return("Number of cycles too low!!!")
  }
  phi<-diag(0/N) ## creating the first phi matrix
  Pt<-t(P)
  x<-0 ## needed for a correct counting cycle
  for (i in 1:n){
    x<-x+1 ## start the counting cycle
    S1<-mtx.exp(S,x) ## powering S
    P1<-mtx.exp(P,x) ## powering P
    Pt1<-mtx.exp(Pt,x) ## powering the transpose of P
    D<-(1-phi)/(2*N) ## calculating the diagonal of the D matrix
    D<-diag(D) ## extracting the diagonal of the above
    D<-diag(D) ## creating the REAL D matix, which is a diagonal matrix
    phi<-phi+(S1%*%Pt1%*%D%*%P1%*%S1) ## Malecot model
  }
  phi
}
