lines.segmented<-function(x, term, bottom=TRUE, shift=TRUE, conf.level=0.95, k=50, 
  pch=18, rev.sgn=FALSE,...){
  if(missing(term)){
          if(length(x$nameUV$Z)>1 ) {stop("please, specify `term'")}
               else {term<-x$nameUV$Z}
               }
  ss<-list(...)
  colore<- if(is.null(ss$col)) 1 else ss$col
  usr <- par("usr")
  h<-(usr[4]-usr[3])/abs(k)
  y<- if(bottom) usr[3]+h else usr[4]-h
  r<- confint.segmented(object=x,parm=term,level=conf.level,rev.sgn=rev.sgn,digits=15)
  m<-r[[term]]
  #FORSE non e' necessaria
  #if(rev.sgn) m<- -m
  #ma invece serve il seguente (se length(psi)=1 e rev.sgn=T):
  m<-matrix(m,ncol=3)
  if(nrow(m)>1) m<-m[order(m[,1]),]
  est.psi<-m[,1]
  lower.psi<-m[,2]
  upper.psi<-m[,3]
  if(length(est.psi)>1) {
      y<- if(shift) y+seq(-h/2,h/2,length=length(est.psi)) else rep(y,length(est.psi))
      }
  segments(lower.psi, y, upper.psi, y, ...)
  points(est.psi,y,type="p",pch=pch,col=colore)
  }
