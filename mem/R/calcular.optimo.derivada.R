calcular.optimo.derivada <-
function(i.curva.map){
  x<-i.curva.map[,1]
  y<-i.curva.map[,2]
  y.s<-loess(y~x)$fitted
  cambio.signo<-abs(diff(sign(diff(diff(y.s)))))
  if (any(cambio.signo!=0)){
    optimo<-1+which.max(cambio.signo)
    resultados<-i.curva.map[x==optimo,]
  }else{
    optimo<-NA
    resultados<-rep(NA,5)
  }
  return(resultados)
}
