t2sGrey<-function(B0,alpha=.38,hct=.4,cbf.base=55,chi=1.8e-7,e.base=.4,w0=267500000,r2=t2Grey(B0,TRUE),relax=TRUE){
  cbv.base <- 0.8 * cbf.base^alpha
  vol.base <- cbv.base/100
  chi.si <- chi * 4 * pi
  w0.si = w0/(2 * pi)
  f.shift.si <- w0.si * B0 * chi.si * hct * (e.base) * 4 * pi/3
  r.si <- vol.base * f.shift.si
  r2.si<-t2Grey(B0,relax=TRUE)
  r2.star.base <- r.si + r2.si
  if(relax){return(r2.star.base)}else{
    return(1/r2.star.base*1000)
  }
}  
