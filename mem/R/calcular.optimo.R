calcular.optimo <-
function(i.curva.map,i.metodo,i.parametro){
  # Metodo 1: Original, segunda derivada
  if (i.metodo==1) resultados<-calcular.optimo.original(i.curva.map,i.parametro)
  # Metodo 2: Usando un criterio de % sobre la pendiente
  if (i.metodo==2) resultados<-calcular.optimo.criterio(i.curva.map,i.parametro)
  # Método 3: Usando la pendiente de la derivada.
  if (i.metodo==3) resultados<-calcular.optimo.pendiente(i.curva.map)
  # Método 4: Segunda derivada, igualando a 0
  if (i.metodo==4) resultados<-calcular.optimo.derivada(i.curva.map)
  # Método 5: Segun Raven
  #if (i.metodo==5) resultados<-calcular.optimo.varadhan(i.curva.map,i.parametro)
  return(resultados)
}
