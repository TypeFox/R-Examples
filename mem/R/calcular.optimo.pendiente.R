calcular.optimo.pendiente <-
function(i.curva.map){
  x<-i.curva.map[,1]
  y<-i.curva.map[,2]
  y.s<-loess(y~x)$fitted
  x.range<-(max(x,na.rm=T)-min(x,na.rm=T))
  y.range<-(max(y,na.rm=T)-min(y,na.rm=T))
  pendiente<-y.range/x.range
  y.d<-diff(y.s)
  optimo<-1+which.min(abs(y.d-pendiente))
  resultados<-i.curva.map[x==optimo,]
  return(resultados)
}
