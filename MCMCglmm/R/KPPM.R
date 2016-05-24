KPPM<-function(m,k){

  if(requireNamespace("combinat", quietly = TRUE)==FALSE){stop("combinat not loaded")}

  K<-matrix(0,m^k, m^k)
  H<-matrix(0,m,m)

  terms<-expand.grid(lapply(1:k, function(x){(1:m)}))
  nterms<-length(terms[,1])

  perms<-combinat::permn(1:k) 
  nperms<-length(perms) 

  for(s.terms in 1:nterms){
    for(k.terms in 1:nperms){
      Ktmp<-1
      for(v.terms in 1:k){ 
        e1<-as.numeric(terms[s.terms,][v.terms])
        e2<-as.numeric(terms[s.terms,][perms[[k.terms]][v.terms]])

        H[,e1][e2]<-1
        Ktmp<-kronecker(Ktmp, H)
        H[,e1][e2]<-0
      }
      K<-K+Ktmp
    }
  }
  K/factorial(k)
}
