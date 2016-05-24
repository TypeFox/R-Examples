kaiserWin<-function(fh=NULL,fl=NULL,tw,sf,d.sa,d.pbr,type){
  
  pbr<-10^(d.pbr/20)-1
  sbr<-10^(d.sa/-20)
  ripple<-min(c(pbr,sbr))
  e.sa<- -20*log10(ripple)
  
  df<-tw/sf
  N<-(e.sa-7.95)/(14.36*df)
  if(round(N)%%2==0){N<-round(N+1)}else{N<-round(N)}
  
  
  if(d.sa>=50)          {B<-.1102*(e.sa-8.7)}
  if(d.sa<50 & d.sa>21) {B<-.5842*(e.sa-21)^.4 + .07886*(e.sa-21)}
  if(d.sa<=21)          {B<-0}
  
  n<- -floor(N/2):floor(N/2)
  k<-B*sqrt(1-(2*n/(N-1))^2)
  
  kaiser.win<-besselI(k,0)/besselI(B,0)
  pulse<-sinc(fh=fh,fl=fl,tw=tw,sf=sf,type=type,n=n)
  
  fir<-pulse*kaiser.win
  return(fir)
}