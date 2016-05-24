memintensity <-
function(i.data,i.levels=c(0.40,0.90,0.975),i.n.max=-1,i.seasons=10){
  datos.modelo<-influenza(i.datos=i.data,i.niveles=i.levels,i.n.max=i.n.max,i.n.max.temp=i.seasons)
  intensity.thresholds<-matrix(c(datos.modelo$pre.post.intervalos[1,3],datos.modelo$epi.intervalos[,4]),ncol=4)
  colnames(intensity.thresholds)<-c("Epidemic",paste(c("Medium (","High (","Very high ("),as.character(datos.modelo$epi.intervalos[,1]*100),"%)",sep=""))
  rownames(intensity.thresholds)<-"Intensity Thresholds"
  results<-list(intensity.thresholds=intensity.thresholds,param.levels=i.levels,param.seasons=i.seasons)
  return(results)
}
