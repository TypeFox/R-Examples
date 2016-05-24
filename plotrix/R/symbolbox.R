symbolbox<-function(x1,y1,x2,y2,tot,relw=0.5,fg=par("fg"),bg=par("bg"),
 box=TRUE,debug=TRUE,...) {

 if(debug) cat("symbolbox:",x1,y1,x2,y2,tot,"\n")
  x <- c(x1,x2)
  y <- c(y1,y2)
  if (x2 < x1) x<-rev(x)
  if (y2 < y1) y<-rev(y)
  pin<-par("pin")
  usr<-par("usr")
  usr.pin<-diff(par("usr"))[c(1,3)]/par("pin")
  dx<-diff(x)/usr.pin[1]
  dy<-diff(y)/usr.pin[2]
  area<-dx*dy
  m<-dx*sqrt(tot/area)
  n<-dy*sqrt(tot/area)
  rm<-max(round(m),1)
  rn<-max(round(n),1)
  while(rm*rn < tot) {
   if((dx*sqrt(tot/area)-m) > (dy*sqrt(tot/area)-n)) {
    rm <- rm + 1
   }
   else {
    rn <- rn + 1
   }
  }
  m<-rm
  n<-rn
  if(debug) cat("symbolbox:",dx,dy,m,n,rm,rn,tot,"\n")
  r<-dx/m*relw/2
  dx<-dx/m*usr.pin[1]
  dy<-dy/n*usr.pin[2]
  mat<-matrix(1:(m*n),nrow=m,ncol=n)
  xpos<-x[1]+(row(mat)[mat <= tot] - 0.5) * dx
  ypos<-y[1]+(col(mat)[mat <= tot] - 0.5) * dy
  symbols(xpos,ypos,rep(1,tot),bg=bg,fg=fg,add=TRUE,inches=r)
  if(box)
   polygon(x[c(1,1,2,2,1)],y[c(1,2,2,1,1)],border=fg,...)
}
