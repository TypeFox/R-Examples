localizar.epidemia <-
function(i.datos,i.n.values=5,i.metodo=2,i.parametro=2.8){
  datos<-as.vector(as.matrix(i.datos))
  curva.map<-calcular.map(datos)
  optimo.map<-calcular.optimo(curva.map,i.metodo,i.parametro)
  if (!is.na(optimo.map[4])){
    pre.epi<-max.n.valores(datos[-(optimo.map[4]:length(datos))],i.n.values)
  }else{
    pre.epi<-max.n.valores(datos,i.n.values)
  }
  if (!is.na(optimo.map[5])){
    post.epi<-max.n.valores(datos[-(1:optimo.map[5])],i.n.values)
  }else{
    post.epi<-max.n.valores(datos,i.n.values)
  }
  return(list(datos=i.datos,curva.map=curva.map,optimo.map=optimo.map,pre.epi=pre.epi,post.epi=post.epi))
}
