#plotmpa: función que grafica las relaciones dentro de un cluster dado
# parámetros: clase: número del cluster que se quiere graficar
#             E: la matriz de asociaciones
#             mpa: la lista resultante de la función mpa
#             fpond: grado de ponderación de los vínculos
#             tam.fuente: tamaño de la fuente
#             tit: título
plotmpa <- function(clase, E, mpa, fpond= 10, tit= NULL, tam.fuente=1)
{
  m <- nrow(E)
  pal.dentro <- seq(1:m)
  npal.dentro <- contar.si(mpa$Clases, clase)
  pal.dentro <- ifelse(mpa$Clases == clase, pal.dentro, NA)
  E1 <- E[!is.na(pal.dentro),!is.na(pal.dentro)]
  if(npal.dentro>1) E1 <- E1-diag(1,nrow(E1),nrow(E1))
  nom <- colnames(E)
  nomE1 <- nom[!is.na(pal.dentro)]
  rownames(E1) <- nomE1
  colnames(E1) <- nomE1
  En <- network(E1)
  color <- rep("blue", times=nrow(E1))
  color1 <- ifelse(nomE1==mpa$Nombres[clase],"red",color)
  plot(En, vertex.cex=2, vertex.sides=100,usearrows=FALSE, displaylabels =TRUE, boxed.label=FALSE, edge.lwd = fpond*E1,
        vertex.col=color1, main=tit, label.cex=tam.fuente)
}
