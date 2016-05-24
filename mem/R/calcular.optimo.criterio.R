calcular.optimo.criterio <-
function(i.curva.map,i.criterio=3){
  y.100<-min((1:length(i.curva.map[,2]))[round(i.curva.map[,2],2)==100])
  curva.map<-i.curva.map[1:y.100,]
  x<-curva.map[,1]
  y<-curva.map[,2]
  y.s<-suavizado(y,1)  
  #y.s<-loess(y~x)$fitted
  d.y<-diff(y.s)
  if (any(d.y<i.criterio)){
    #optimo<-max((1:(length(x)-1))[d.y>i.criterio],na.rm=T)
    optimo<-min((1:(length(x)-1))[d.y<i.criterio],na.rm=T) 
  }else{
    optimo<-length(d.y)+1
  }
  resultados<-curva.map[x==optimo,]
  return(resultados)
}
