#Funcion auxiliar para grafica bivariada con marginales

trans3d<-function (x, y, z, pmat) 
{
    tr <- cbind(x, y, z, 1) %*% pmat
    list(x = tr[, 1]/tr[, 4], y = tr[, 2]/tr[, 4])
}