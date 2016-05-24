prev<-function(f0,f1,f2,g,tmother,tfather,maf){
  dec=3
  tm=tmother
  tf=tfather
  p=maf

  mu1=(p^4)
  mu2=2*(p^3*(1-p))
  mu3=(p^2*(1-p)^2)
  mu4=4*(p^2*(1-p)^2)
  mu5=2*(p*(1-p)^3)
  mu6=((1-p)^4)

  tau=c(1,
       tf,(1-tf),tm,(1-tm),
       1,1,
       tm*tf,
       tm*(1-tf),(1-tm)*tf,
       (1-tm)*(1-tf),
       tm,(1-tm),tf,(1-tf),
       1)

  f=c(f2,
      f2,f1*g,f2,f1,
      f1*g,f1,
      f2,f1*g,f1,f0,
      f1*g,f0,f1,f0,
      f0)

  mu<-c(mu1,rep(mu2,4),rep(mu3,2),rep(mu4,4),rep(mu5,4),mu6)

  d=signif(sum(f*mu*tau),dec)
  names(d)='d'
  d
}
