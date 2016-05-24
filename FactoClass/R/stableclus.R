stableclus= function(dudi,part,k.clust,ff.clus=NULL,bplot=TRUE,kmns=FALSE) 
{
    if(k.clust<2 || part<2) {stop("Number of partitions #Condicionales a los requisitos de la funcion
        and classes must be more than 2")
    }
    if((k.clust^part)>100000) {stop("Please select 
        a smaller number of partitions or clusters")
    } 
    if(class(dudi)[2]!="dudi"){stop("Not a valid dudi object")
    } 

    nf=dudi$nf                    #Extracción de los objetos necesarios del
    obj.clasf=dudi$li             #objeto dudi
    pesos=dudi$lw
    n=nrow(obj.clasf)  
    A=matrix(c(rep(0,part*(k.clust^part))),ncol=part)         
    m=nrow(A)
    cont=(k.clust^(part-1))
    i=ncol(A)-1
    c=c(1:k.clust)
    k=k.clust
    A[,part]=rep(c,cont)   #build up the matrix 
    while(cont>1){          #to save information 
        t=c(rep(k.clust,k)) #of the all classifications                                    
        c=rep(c,t)          #of the one individual
        cont=cont/k.clust
        A[,i]=rep(c,cont)
        i=i-1
        k=k.clust*k
    }
    ID=c(1:m)
    cluster=matrix(c(rep(0,n*part)),ncol=part)                 #Aqui se guarda la informacion
    for (i in 1:part) {                                        #de la clasificacion, en cada una
        kmeans=as.vector(kmeansW(x=obj.clasf,centers=k.clust,  #de las particiones, de todos los 
        weight=pesos)$cluster)                                 #individuos
        cluster[,i]=kmeans  
    }
    f=c(rep(0,m))              #En esta parte se cuentan el numero de individuos que pertenecen a
    ide=c(rep(0,n))            #cada una de las nuevas clases producto y se guardan estas frecuencias
    for (i in 1:m){
        for(j in 1:n){
            if(identical(cluster[j,],A[i,]) ){ 
                f[i]=f[i] + 1
                ide[j]=i
            } 
        }
    }  
    ide2=c(1:n)
    l=order(f,decreasing=TRUE) #Se procede a ordenar las frecuencias y     
    fo=sort(f,decreasing=TRUE) #se presentan en un diagrama de barras 
    fot=fo[fo>0]               #para que el usuario decida el número de clases
    IDo=ID[l]                  #finales
    IDot=IDo[fo>0]
    if(bplot){
        barplot(fot)
    }
    if(ff.clus=="NULL"){  
        cat("Select the number of clusters:")
        ff.clus <- as.integer(readLines(n=1))
    }
    IDotf=IDot[1:ff.clus]
    IDotff=factor(IDotf)
    clsfrts=length(IDotff)                   #Ahora se procede a calcular los centros de gravedad de las                     
    cls.inc=list()                           #clases seleccionadas
    for(i in 1:clsfrts){
        cls.inc[[i]]=as.matrix(obj.clasf[ide==IDotff[i],])
    }
    f1=function(X){ tapply(X,col(X),mean) }
    c.grav=lapply(cls.inc,f1)
    C.grav=matrix(0,ncol=nf,nrow=ff.clus) 
    for(i in 1:ff.clus){
        C.grav[i,]=c.grav[[i]]
    }
    val=c(rep(0,n))
    for (i in 1:ff.clus){
        for (j in 1:n){        
            if(ide[j] != IDotff[i]){val[j]=val[j]+1}
        }   
    }         
    Reafct=obj.clasf[val==ff.clus,]            #Finalmente se agregan los demas individuos
    ide3=ide2[val==ff.clus]                    #por reafectacion a las clases definitivas
    n.reafct=nrow(Reafct)                      #para esto se calculan las distancias a cada uno
    fdist=matrix(0,ncol=ff.clus,nrow=n.reafct) #de los centros de gravedad para cada individuo y se clasifica según
    for (i in 1:ff.clus){                      #la distancia minima
        dist=(Reafct-matrix(C.grav[i,],nrow=n.reafct,ncol=nf,byrow=TRUE))^2
        fdist[,i]=as.vector(sqrt(tapply(dist,row(dist),sum)))
    } 
    class=c()
    for(i in 1:n.reafct){  
        for(j in 1:ff.clus){
            if(fdist[i,j]==min(fdist[i,])){class[i]=IDotf[j]}
        }
    }
    ide[ide3]=class                           
    ide=as.factor(ide)
    ide=as.numeric(ide)
    ide 
    if(kmns){
        ide=kmeansW(x=obj.clasf,centers=C.grav)$cluster
        C.grav=kmeansW(x=obj.clasf,centers=C.grav)$centers
    } 
    
    OUTPUT <- list(cluster = ide, centers = C.grav ) #El resultado es una lista con
    class(OUTPUT) <- "stableclus"                    #los centros de gravedas y el vector
    return(OUTPUT)                                   #que identifica las clases
                                          
                                         
                                         
}
  