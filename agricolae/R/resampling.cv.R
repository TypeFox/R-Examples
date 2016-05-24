`resampling.cv` <-
function(A,size,npoints) {
nc.A<-ncol(A)
orden<-1:nc.A
y<-trunc(seq(2,(nc.A-1),length=npoints))
marcador<- unique(y)
nmarca<-length(marcador)
CV<-rep(0,nmarca)
# Control del tiempo inicio y final
inicio<-Sys.time()
i<-0
# Proceso de re-muestreo
for (k in marcador) {
cv.m <- 0
for (m in 1:size) {
muestra <-sample(orden,k)
B<-A[,muestra]
cv.m <- cv.m + cv.similarity(B)
}
i<-i+1
CV[i]<- cv.m/size
}
# Final del proceso de re-muestreo
final<-Sys.time()
# Impresion de resultados
time<-final-inicio
cat("\nTime of process ...",time,"\n")
tabla.cv <- data.frame(marcador=marcador,CV)
# Estimacion de un modelo para CV.
modelo<-lm(CV ~ I(1/marcador))
return(list(model=modelo,table.cv=tabla.cv))
}

