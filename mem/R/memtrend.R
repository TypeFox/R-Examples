memtrend <-
function(i.data,i.seasons=10){
  temporadas<-dim(i.data)[2]
  if (i.seasons<1) i.seasons<-temporadas
  temporadas.usar<-(max(temporadas-i.seasons+1,1)):temporadas
  n.temp<-length(temporadas.usar)
  datos<-i.data[temporadas.usar]
  datos.2<-apply(datos,2,fill.missing)
  datos.3<-apply(datos.2,2,diff)
  datos.3<-rbind(NA,datos.3)
  rownames(datos.3)[1]<-rownames(datos.2)[1]
  datos.flu<-influenza(datos,i.n.max.temp=-1)
  for (j in 1:ncol(datos.3)){
    i.i<-datos.flu$gripe[1,j,1]-1
    i.f<-datos.flu$gripe[2,j,1]+2
    if (i.i>=1) datos.3[1:i.i,j]<-NA
    if (i.f<=nrow(datos.3)) datos.3[i.f:nrow(datos.3),j]<-NA
  }  
  datos.4<-datos.3
  datos.5<-datos.3
  datos.4[datos.3<0]<-NA
  datos.5[datos.3>0]<-NA
  limite.s<-iconfianza.completo(as.numeric(datos.4),tipo=1)$mayor[1]
  limite.i<-iconfianza.completo(as.numeric(datos.5),tipo=1)$menor[2]
  trend.thresholds<-matrix(c(limite.s,limite.i),ncol=2)
  colnames(trend.thresholds)<-c("Ascending Threshold","Descending Threshold")
  rownames(trend.thresholds)<-"Trend Thresholds"
  results<-list(trend.thresholds=trend.thresholds,param.seasons=i.seasons)
  return(results)
}
