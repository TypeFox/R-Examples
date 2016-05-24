calcular.optimo.original <-
function(i.curva.map,i.hsuav){
  x<-i.curva.map[,1]
  y<-i.curva.map[,2]
  y.d<-diff(y)
  x.d<-x[-length(x)]
  #y.s<-loess(y.d~x.d)$fitted
  y.s<-suavizado(y.d,0.6)
  x.n<-normalizar(x.d)
  y.n<-normalizar(y.s)
	u<-(x.n-y.n)/sqrt(2)
	v<-sqrt(x.n^2+y.n^2-u^2)
  optimo<-which.min(v)
  resultados<-i.curva.map[x==optimo,]
  return(resultados)
}
