#Calculation of simultaneous confidence intervals for bifactorial
#designs based on the multivariate t-distribution
intstudent2<-function(C,nboot,simerror,alpha,...){
  dauer<-proc.time()[3]
  groups<-carpetkombi(C)
  y<-carpetdaten(C)
  contrasts<-carpetmatrix(C)
  nctr<-length(contrasts[,1])
  intervals<-confint(glht(lm(y~groups),linfct=contrasts))$confint
  kiu<-intervals[nctr+1:(2*nctr)]
  kio<-intervals[(2*nctr)+1:(3*nctr)]
  cnames<-character(0)
  for(a in 1:C@D[1]){for(b in 1:C@D[2]){
    cnames<-c(cnames,paste("(",a,",",b,")-(",a,",0)",sep=""),paste("(",a,",",b,")-(0,",b,")",sep=""))
  }}
  dauer=proc.time()[3]-dauer
  new("margint",cnames=cnames,kiu=kiu,kio=kio,method="Multivariate t-distribution",nboot=0,
      simerror=0,alpha=alpha,duration=dauer)
}
