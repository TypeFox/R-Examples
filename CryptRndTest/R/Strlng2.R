Strlng2=function(n,k,log=TRUE){
  # Kaynaklar: 
  # W W Bleick and Peter C. C. Wang, Asymptotics of Stirling Numbers of the Second Kind, 
  #        Proceedings of the American Mathematical Society Vol. 42, No. 2 (Feb., 1974), pp. 575-580
  # N M Temme, Asymptotic Estimates of Stirling Numbers, STUDIES IN APPLIED MATHEMATICS,
  #        89:233-243 (1993), Elsevier Science Publicshing.
  # This function requires the R package: LamberW.
  if (k==1){
    topl=1
  }else{
    z=1:n
    nu=n/k
    r=1:(n-k)
    G=-W(-nu*exp(-nu))
    topl=log(sqrt(n-k))+(n-k)*log((n-k)/exp(1))+
        (sum(log(z))-sum(log(r))-sum(log(1:k)))-log(sqrt(n*(1-G)))-k*log(G)-(n-k)*log(nu-G)
  }
  if (log==FALSE){
     topl=exp(Rmpfr::mpfr(topl,250))     
  }
  return(Stirling.num=topl)
}
