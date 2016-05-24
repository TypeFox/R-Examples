#diagram.mpa: función que grafica el diagrama estratégico de los grupos
#             mpa: la lista resultante de la función mpa
#             tit: título
diagram.mpa <- function(mpa, tmin=3, tit= NULL, pos=1)
{
  x <- mpa$Resumen[,4]
  y <- mpa$Resumen[,3]
  nom <- mpa$Nombres
  tam <- mpa$Resumen[,2]
  x <- x[tam>=tmin]
  y <- y[tam>=tmin]
  nom <- nom[tam>=tmin]
  plot(x, y, type="p", main=tit, xlab="Centralidad", ylab="Densidad",las=1)
  abline(h=mean(y, na.rm=TRUE),col="grey");  abline(v=mean(x, na.rm=TRUE),col="grey")
  text(x, y, labels=nom, pos=pos)
}
