fSombra6<-function(angGen, distances, struct, prom=TRUE)
{
    stopifnot(is.list(struct),
              is.data.frame(distances))
    ##distances sólo tiene tres distances, así que genero una cuadrícula
    if (dim(distances)[1]==1){ 
        Red<-with(distances,
                  data.frame(Lew=c(-Lew,0,Lew,-Lew,Lew),
                             Lns=c(Lns,Lns,Lns,0,0),
                             H=H))
    } else { #distances es una matriz, luego no hace falta generar la cuadrícula
        Red<-distances[1:5,]} #Sólo necesito las 5 primeras filas...necesario por sí  se entrega un data.frame erroneo

    ## Calculo la sombra debida a cada uno de los 5 seguidores
    SombraGrupo<-matrix(ncol=5,nrow=dim(angGen)[1])
    for (i in 1:5) {SombraGrupo[,i]<-coredata(fSombra2X(angGen,Red[i,],struct))}
    ##Para calcular la Sombra Promedio, necesito el número de seguidores en cada posición(distrib)
    distrib=with(struct,c(1,Ncol-2,1,Nrow-1,(Ncol-2)*(Nrow-1),Nrow-1)) 
    vProm=c(sum(distrib[c(5,6)]),
            sum(distrib[c(4,5,6)]),
            sum(distrib[c(4,5)]),
            sum(distrib[c(2,3,5,6)]),
            sum(distrib[c(1,2,4,5)]))
    Nseg=sum(distrib) ##Número total de seguidores
    ##Con la función SWEEP multiplico el Factor de Sombra de cada tipo (columnas de SombraGrupo) por el resultado de vProm
    
    if (prom==TRUE){
        ## Factor de Sombra Promedio en el grupo de SEIS seguidores teniendo en cuenta distrib
        FS=rowSums(sweep(SombraGrupo,2,vProm,'*'))/Nseg
        FS[FS>1]<-1
    } else {		
        ## Factor de sombra en el seguidor #5 debido a los otros 5 seguidores
        FS=rowSums(SombraGrupo)
        FS[FS>1]<-1}
    return(zoo(FS, index(angGen)))
}
