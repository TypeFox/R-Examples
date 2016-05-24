lk_ar_rho <- function(lrho,SUP,V,outp=FALSE){
  q = ncol(SUP)
  rho = expit(lrho)*2-1
  si2 = 1-rho^2
  lWei = -(t(SUP)-rho*SUP)^2/si2/2-log(2*pi*si2)/2
  lWei = lWei-log(rowSums(exp(lWei)))%o%rep(1,q)
#  Wei = normpdf(SUP',rho*SUP,sqrt(1-rho^2));
#  Wei = diagv(1./sum(Wei,2),Wei);
#  Wei1 = max(Wei,10^-100);
#  lk = sum(sum(V.*log(Wei1)));
  lk = sum(V*lWei)  
  flk = -lk
  if(outp) flk = list(flk=flk,Wei = exp(lWei),rho=rho)
  flk
  
}