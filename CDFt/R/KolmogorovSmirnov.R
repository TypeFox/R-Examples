KolmogorovSmirnov <- function(S1, S2){

  xS1=sort(S1)
  cdftmp=ecdf(xS1)
  cdf1=cdftmp(xS1)

  xS2=sort(S2)
  cdftmp=ecdf(xS2)
  cdfEstim=cdftmp(xS2)
      
  cdfRef=approx(xS1,cdf1,xS2,yleft=0,yright=1, ties='mean')
  
  dif=cdfRef$y-cdfEstim
  dif=abs(dif)
      
  Ks=max(dif)
  
  return(Ks)

}
