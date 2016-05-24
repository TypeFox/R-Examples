#############################################################################################
##  Funcion de enlace: Combinación de métodos factoriales y de clasificaciòn no            ##
##  supervisada.                                                                           ##
##                                                                                         ##
##  Mayo 15 de 2011 inclusión de parámetro de pesos en función FactoClass (CEPT)           ##
##                                                                                         ##
## Elaborado por: Pedro Cesar del Campo Neira                                              ##
## Revisado y modificado por: Campo Elías Pardo  INGLÉS  Nov.30/07                         ##
## Universidad Nacional de Colombia                                                        ##
##                                                                                         ##
## requiere:ade4      library(ade4)                                                        ##
##                                                                                         ##
## Fac.Num  ( dfact   = objeto 'data.frame' de datos variables activas,                    ##
##            metodo  = funcion de ade4 para metodo factorial.                             ##
##            dfilu   = variables ilustrativas (deafault NULL)                             ##
##            nfaf    = Numero de ejes para el analisis (deafault 2)                       ##
##            nfcl    = Numero de ejes para la clasificación (deafault NULL)               ##
##            k.clust = Numero de clases (deafault NULL)                                   ##
##            scanFC  = 'TRUE',escanea ,y si 'FALSE', no escanea                           ##
##            n.max   = si 'dim(dfact)[1]>=n.max' efectua previo k-means (deafault 5000)   ##
##            n.clus  = si 'dim(dfact)[1]>=n.max' efectua WARD con n.clus (deafault 1000)  ##
##            sign    = valor estadistico de rechazo en las pruebas.                       ##
##            conso   = realiza proceso de consolidación de la clasificación(deafault TRUE)##
##            n.indi  = número de indices en el grafico (default 25)                       ##
##          )                                                                              ##
##                                                                                         ##
#############################################################################################


FactoClass<-function( dfact, metodo,dfilu = NULL , nf = 2, nfcl = 10, k.clust = 3, 
                      scanFC = TRUE , n.max = 5000 , n.clus = 1000 ,sign = 2.0,
                      conso=TRUE , n.indi = 25,row.w = rep(1, nrow(dfact)))
{

  n <- dim(dfact)[1]

  n.act  <- deparse(substitute(dfact))  ### Tipo caracter nombre de dfact
  metodo <- deparse(substitute(metodo)) ### Tipo caracter nombre de la función
 ### construccion del llamado función dudi
  row.w <- row.w/sum(row.w) # asegurar que los pesos suman 1
  if(metodo=="dudi.coa") call1 <- call(metodo,df = as.name(n.act), nf = nf , scannf = scanFC)
     else call1 <- call(metodo,df = as.name(n.act), nf = nf , scannf = scanFC,row.w=row.w) 
                                                                        
  par(las=1)                                                                      
  DuDi1 <- eval(call1) # evaluación del llamado función dudi.*
  nf    <- DuDi1$nf
  cat("The number of retained axes for factorial analysis is ",nf,"\n\n") 

  if(scanFC==TRUE){  #### Selecciona numero de ejes para realizar el proceso de clasificación
    cat("Select the number of axes for clustering: ")
    nfcl <- as.integer(readLines(n = 1))
  }

  DuDi2 <- redo.dudi( DuDi1, newnf = nfcl ) ### objeto dudi para clasificación
  nfcl <- DuDi2$nf
 
  cat("The number of axes for clustering is ",nfcl,"\n\n")
             
  objetos   <- DuDi2$li  ### ejes factoriales de filas para clasificación
  pesos     <- DuDi2$lw  ### pesos de filas para clasificación
  obj.clasf <- objetos   ### elementos que entran a la clasificación

###########################################################################
######################### Primer criterio de clasificacion "n >= n.max"

  if(n >= n.max){      	
    prev.kmeans <- kmeansW(x = obj.clasf, centers = n.clus, weight = pesos)
    obj.clasf   <- prev.kmeans$centers
    pesos       <- tapply(pesos, prev.kmeans$cluster, sum)
    prev.size   <- prev.kmeans$size
  }

###########################################################################
######################### clasificación no supervisada método de WARD 

  dend <- ward.cluster( dista= dist(obj.clasf), peso=pesos ,h.clust = 0, n.indi = n.indi)
  cat("Look the histogram of",n.indi,"indexes \n")

  if(scanFC == TRUE){### Selecciona numero el número de clases
    cat("Select the number of clusters: ")       
    k.clust <- as.integer(readLines(n = 1))
  }
  cat("Partition in ", k.clust, " clusters\n")

  cluster1 <- cutree(dend$HW, k = k.clust) ### Clasificacion generada por WARD

  if(n >= n.max){
    d1 <- data.frame(prev = prev.kmeans$cluster, id = 1:n)
    d2 <- data.frame(prev = names(cluster1), cl2 = cutree(dend$HW, k = k.clust))
    
    dd <- merge(d1, d2, all.x = TRUE)
    dd <- dd[order(dd$id), ]
    cluster1 <- dd$cl2
  }

  dev.new()     ### Dendograma con clasificación
  plot(dend$HW,las=1,sub="",xlab="",ylab="Indexes",main="")
  rect.hclust(dend$HW, k.clust, border="blue")

###########################################################################
######################### K-MEANS

  ft  <- function(x){data.frame(t(x))}  # funcion transpone y convierte en data.frame
  pes <- function(x){x/sum(x)}          # funcion para convertir en peso de cada clase

###----------------------------------------------------------------------------------

  p.clust1 <- DuDi2$lw
  for (k in 1:k.clust){ p.clust1[cluster1==k] <- pes(p.clust1[cluster1==k]) } ## Pesos de los individuos 
                                                                              ## para cluster 1.                                                                            
  center1 <- lapply(by( p.clust1 * DuDi2$li , cluster1, colSums ),ft)         ## Centros de la clasificación 
  center1 <- list.to.data(center1)[-1]                                        ## generada por WARD 

#################### ordena la clasificación por el primer componente principal

#  critrio.orden     <- order(center1[,1])
#  center1           <- center1[ critrio.orden, ]
#  rownames(center1) <- NULL
#  
#  d1 <- data.frame(prev = cluster1  , id = 1:n)
#  d2 <- data.frame(prev = 1:k.clust , ordenado = critrio.orden)    
#  dd <- merge(d1, d2, all.x = TRUE)
#  dd <- dd[order(dd$id), ]
#  cluster1 <- dd$ordenado
  
  
 
 if(conso){    ########################################### con consolidación   
    clus.summ <- NULL 
###########################################################################
###########################################################################

  ###  clasificación generada por K-MEANS con centros de WARD(center1)
    cluster2 <- kmeansW( x = objetos , centers = center1 , weight = pesos )$cluster    
    #cluster2 <- kmeans( objetos , center1)$cluster 
                   
###########################################################################
######################### PROPIEDADES DE LA CLASIFICACION (cluster2)
                                             ## Tabla de comportamiento de inercia de las clases 2
  for(k in 1:k.clust ){
    clus.summ <- rbind( clus.summ , analisis.clus(DuDi2$li[cluster2==k,],DuDi2$lw[cluster2==k]) ) 
  }
  
  clus.summ <- rbind( clus.summ , apply(clus.summ,2,sum)  )
  clus.summ[k.clust + 1,4] = NA
  rownames(clus.summ)[k.clust + 1] <- "TOTAL"

  clus.summ1 <- NULL

    for( k in 1:k.clust ){
      clus.summ1 <- rbind( clus.summ1 , analisis.clus(DuDi2$li[cluster1==k,],DuDi2$lw[cluster1==k]) )  
    }
      
      clus.summ1 <- rbind( clus.summ1 , apply(clus.summ1,2,sum)  )
      clus.summ1[k.clust + 1,4] <- NA  
      clus.summ <- data.frame( Bef.Size       =   clus.summ1$Size     ,
                               Aft.Size       =   clus.summ$Size      ,
                               Bef.Inertia    =   clus.summ1$Inertia  ,
                               Aft.Inertia    =   clus.summ$Inertia   ,
                               Bef.Weight     =   clus.summ1$Weight   ,
                               Aft.Weight     =   clus.summ$Weight    ,
                               Bef.Dist_2     =   clus.summ1$Dist_2   ,
                               Aft.Dist_2     =   clus.summ$Dist_2)

  rownames(clus.summ)[k.clust + 1] <- "TOTAL"  
    
  } # fin consolidación                                                 

  if(!conso){########### --------------------- sin consolidación
    clus.summ1 <- NULL
               ## Tabla de comportamiento de inercia de las clases 1 y 2
      for(k in 1:k.clust){
        clus.summ1 <- rbind(clus.summ1 , analisis.clus(DuDi2$li[cluster1 == k, ], DuDi2$lw[cluster1 == k]))
      }
                                      
      clus.summ <- data.frame( Size             =   clus.summ1$Size     ,
                               Inertia          =   clus.summ1$Inertia  ,
                               Weight           =   clus.summ1$Weight   ,
                               Dist_2           =   clus.summ1$Dist_2)
                         
      clus.summ <- rbind( clus.summ1 , c( sum(clus.summ[1]) ,
                                          sum(clus.summ[2]) ,
                                          sum(clus.summ[3]) ,
                                          NA))
      rownames(clus.summ)[k.clust + 1] <- "TOTAL"
    
  }
###########################################################################
######################### COORDENADAS DE LAS CLASES (cluster2)
  if (!conso) cluster2 <- cluster1             
    
  p.clust <- DuDi1$lw 
  for (k in 1:k.clust)p.clust[cluster2==k] <- pes(p.clust[cluster2==k]) 

  cor.clus <- lapply(by( p.clust * DuDi1$li , cluster2, colSums ),ft)
  cor.clus <- list.to.data(cor.clus)[-1]

###########################################################################
######################### CARACTERIZACION DE LA CLASIFICACION (cluster2)

  base0 <- dfact
  
#if(class(DuDi1)[1] == "coa" ){ base0 <- data.frame(t(t(dfact)/colSums(dfact))) }

  if( is.null(dfilu) == FALSE ){ 
   if(class(dfilu)!="data.frame"){ return(cat("\n\n ERROR: Illustrative Variables should be 'data.frame'\n")) }
   if(dim(dfilu)[1]!= n ){ return(cat("\n\n ERROR: Active and  Illustrative Variables 
                           should have the same number of elements\n")) }
   base0 <- data.frame(base0,dfilu) 
  }

  base0 <- Fac.Num(base0)

  carac.cont = NULL
  carac.cate = NULL
  carac.frec = NULL
  carac.fril = NULL

  if(is.null(base0$numeric)==FALSE){ carac.cont <- cluster.carac( base0$numeric, cluster2 ,"co", sign) }
  if(is.null(base0$factor )==FALSE){ carac.cate <- cluster.carac( base0$factor , cluster2 ,"ca", sign) } 
# agregado por CEPT mayo 14/09
  if(is.null(base0$integer)==FALSE){ carac.frec <- cluster.carac( base0$integer , cluster2 ,"fr", sign) }
  
  if(class(DuDi1)[1] == "coa" ){
    if(is.null(dfilu)==FALSE) dfact <- data.frame(dfact,dfilu)
    carac.frec <- cluster.carac(dfact,cluster2,"fr",sign)
  }
###########################################################################
###########################################################################

  cluster2 <- factor(cluster2)

###########################################################################
######################### SALIDA 

  SALIDA <- list( dudi2      = DuDi2, 
                  dudi       = DuDi1,
                  nfcl       = nfcl,
                  k          = k.clust,
                  indices    = dend$INDICES,
                  cluster    = cluster2,
                  cor.clus   = cor.clus,
                  clus.summ  = clus.summ,
                  carac.cont = carac.cont,
                  carac.cate = carac.cate,
                  carac.frec = carac.frec )

  class(SALIDA) <- "FactoClass"

  return(SALIDA)

}
####################################################################################################
#########################          FIN DEL PROGRAMA        #########################################
####################################################################################################
























####################################################################################################
#########################     FUNCION DE ANALISIS EN CLUSTER     ###################################
####################################################################################################

analisis.clus <- function(X,W){

 si <- dim(X)[1]
 Wo <- round( W/sum(W)               , 4 )
 mX <- colSums(Wo*X)
 Xc <- t(t(X)-mX)
 We <- round( sum(W)                 , 4 )
 di <- round( sum(mX^2)              , 4 )

 In <- round( sum(rowSums(Xc^2)* W ) , 4 )
 
 SALIDA <- data.frame(Size     = si, 
                      Inertia  = In,
                      Weight   = We,
                      Dist_2   = di
                     )
 
 return(SALIDA)

}
####################################################################################################


