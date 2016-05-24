####################################################################################################################
####################################################################################################################
####################################################################################################################
###                                                                                                              ###
### CLASIFICACION JERARQUICA METODO DE WARD                                                                      ###
###                                                                                                              ###
### ## Elaborado por: Pedro Cesar del Campo Neira                                                                ###
###    Revisado por: Campo Elías Pardo                                                                           ###
###    Universidad Nacional de Colombia                                                                          ###
###                                                                                                              ###
### ward.cluster( dista  := matriz de distancias euclidianas de los elementos                                    ###
###               peso   := vector de pesos de los elementos (Asume pesos iguales)                               ###
###               plots  := (TRUE O FALSE) muestra un dendograma e indices                                       ###
###               h.clust:= '0' devuelve salida tipo 'hclust' y tabla de Indices                                 ###
###                         '1' devuelve salida tipo 'hclust'                                                    ###
###                         '2' devuelve tabla de Indices                                                        ###
###               n.indi := número de índices a grficar en el histograma (default 25)                            ###                 ###      
###             )                                                                                                ###
###                                                                                                              ###
####################################################################################################################
####################################################################################################################

ward.cluster <- function(dista, peso = NULL , plots = TRUE, h.clust = 2, n.indi = 25 ){


n <- as.integer(attr(dista, "Size"))        # Cantidad de elementos dados por dista
distaM <- as.matrix(dista)                  # dista como matriz


if(is.null(peso)==TRUE){ peso <- rep(1,n) } # Pesos iguales cuando (peso = NULL)

peso=peso/sum(peso)                         # ponderacion de suma 1


fw <-function(a,b){(a*b)/(a+b)}             # funcion ponderación pesos inicial de Ward

distW <- distaM^2 * outer(peso,peso,"fw")   # Matriz inicial en metodo de Ward
distW <- as.dist(distW)             # Matriz inicial en metodo de Ward tipo dist


HW    <- hclust(distW, method="ward.D", members=peso)

#-------------------

   if(h.clust==1){return(HW)}

#-------------------


if(plots==TRUE){                      #Grafico dendograma e histograma
                
                dev.new()
                par(las=1)
                plot(HW,las=1,sub="",xlab="",ylab="Indexes",main="")
                
                dev.new()
                histog <- HW$height[order(HW$height,decreasing=TRUE)]
                histog <- histog[1:n.indi]
                par(las=1)
                barplot(histog,horiz=TRUE)
               }

#-------------------
  
  Nodo <- ( 1:(n-1) ) + n    # Nodo
  Prim <- HW$merge[,1]       # Primogenito
  Benj <- HW$merge[,2]       # Benjamín
  
  
  SALIDA <- data.frame(Nodo,Prim,Benj)
  
  SALIDA[SALIDA[,2]>0,2] <- SALIDA[SALIDA[,2]>0,2] + n 
  SALIDA[SALIDA[,2]<0,2] <- abs(SALIDA[SALIDA[,2]<0,2] )  # Arreglo  Primogenito
  
  SALIDA[SALIDA[,3]>0,3] <- SALIDA[SALIDA[,3]>0,3] + n
  SALIDA[SALIDA[,3]<0,3] <- abs(SALIDA[SALIDA[,3]<0,3] )  # Arreglo  Benjamín
  
  SALIDA[,1] <- factor(SALIDA[,1])
  SALIDA[,2] <- factor(SALIDA[,2])                        # Arreglo a factores
  SALIDA[,3] <- factor(SALIDA[,3])
  
  SALIDA <- data.frame(SALIDA, Indice = HW$height )       # Agregando indice
  
  
  if(h.clust==0){return(list(HW=HW,INDICES=SALIDA))}
  if(h.clust==2){return(SALIDA)}

}
####################################################################################################
############ FIN DE LA FUNCION #####################################################################
####################################################################################################################
####################################################################################################################


