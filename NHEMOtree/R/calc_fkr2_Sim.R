calc_fkr2_Sim <-
function(var_index, ZHT, Varis, CO, index_set, Daten, cut){
  queue          <- list()
  queue[[1]]     <- list()
  queue[[1]][[1]]<- var_index
  queue[[1]][[2]]<- index_set
  queue[[1]][[3]]<- cut
  fk             <- 0
  
  while(length(queue)>0){
    vari_index<- queue[[1]][[1]]

    # Indizes der an diesem Knoten aufzuteilenden Daten (Beobachtungen) 
      akt_index_set<- queue[[1]][[2]]    
      cut          <- queue[[1]][[3]]
      queue[[1]]   <- NULL
    
    # Linke und rechte Knoten bestimmen
      wurzel<- ZHT[which(ZHT[,1]==vari_index),]
      left  <- wurzel[which(wurzel[,2]==0),]
      right <- wurzel[which(wurzel[,2]==1),]
    
      variable<- Varis[[vari_index,2]]
    
    # Linke und Rechte Teilmenge der Beobachtungen berechnen
      left_index_set <- intersect(akt_index_set, which(Daten[, variable] <= cut)) 
      right_index_set<- intersect(akt_index_set, which(Daten[, variable] >  cut))
    
    # Fuege linken Knoten ggf. der Queue hinzu
      if(left[3]!=-99 & length(left_index_set)!=0){
        queue[[(length(queue)+1)]]<- list(left[3], left_index_set, CO[left[3],2])
      }else if(left[3]==-99){ # Wenn Blatt berechne Fehlklassifikation
        class   <- Daten[left_index_set,1]
        fk_local<- length(class)-max(table(class))
        fk      <- fk+fk_local
      }
      
    # Fuege rechten Knoten ggf. der Queue hinzu        
      if(right[3]!=-99 & length(right_index_set)!=0){
        queue[[length(queue)+1]]<- list(right[3], right_index_set, CO[right[3],2])
      }else if(right[3]==-99){
        class   <- Daten[right_index_set,1]         # Klasse steht im Simulationsdatensatz immer in 1. Spalte
        fk_local<- length(class)-max(table(class))
        fk      <- fk+fk_local
      }
  } # Ende der while-Schleife

  # Ausgabe:
    return(fk)
}
