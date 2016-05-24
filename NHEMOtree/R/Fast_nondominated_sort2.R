Fast_nondominated_sort2 <-
function(FKR, Kosten, IDs, popsize){

  Fitness_Original<- rbind(FKR, Kosten)
  i               <- order(Fitness_Original[1,], Fitness_Original[2,])
  Fitness_Original<- Fitness_Original[,i]  

  # Sortieren nach Raengen fuer Original-Baum
    temp_O                 <- rbind(nds_rank(Fitness_Original), IDs, Fitness_Original)
    ii                     <- order(temp_O[1,], temp_O[3,], temp_O[4,], temp_O[2,])
    Fitness_sorted_Original<- temp_O[,ii]           

  # Rang an Position All/2
    cd_rank      <- Fitness_sorted_Original[1,popsize]
    which_cd_rank<- which(Fitness_sorted_Original[1,]<=cd_rank)
    n_toomany    <- length(which_cd_rank)-popsize

  # Falls Crowding Distance Sorting NICHT notwendig
    if (n_toomany == 0){
        Fitness_sorted_Original<- Fitness_sorted_Original[,which(Fitness_sorted_Original[1,]<=cd_rank)]
    }

  # Falls Crowding Distance Sorting notwendig
    if (n_toomany > 0){
      # Alle Individuen mit Rang <= 'Temp'-Rang
        FS_after_NDS<- Fitness_sorted_Original[,which_cd_rank]

      # Berechnung der Crowding-Distance fuer alle Individuen mit Rang 'cd_rank'
        Crowding_Set<- FS_after_NDS[,which(FS_after_NDS[1,]==cd_rank)]
        n_keep      <- ncol(Crowding_Set)-n_toomany

      # Falls Crowding-Set nur aus zwei Individuen besteht, waehle zufaelliges eines der beiden aus!
        if (ncol(Crowding_Set)==2){
	      temp<- runif(1,0,1)
	      if (temp <  0.5) Crowding_Set_kept<- Crowding_Set[,1]
	      if (temp >= 0.5) Crowding_Set_kept<- Crowding_Set[,2]
        }

      if (ncol(Crowding_Set) > 2){
      # Berechne Distanzmasse fuer die m Zielfunktionen
        Distance_overall<- matrix(rep(0, ncol(Crowding_Set)), ncol=1)

        for (m in 1:2){
        # Sortieren
          ii              <- order(Crowding_Set[(m+2),])
          Crowding_Set    <- Crowding_Set[,ii]
          Distance_overall<- t(as.matrix(Distance_overall))[,ii]
          Distance        <- rep(0, ncol(Crowding_Set))
          Crowding_Set    <- rbind(Crowding_Set, Distance)

        # Berechne Distanzmass bzgl. m. Zielfunktion fuer alle Werte, die keine Extremwerte sind.
          Spannweite<- Crowding_Set[(m+2), ncol(Crowding_Set)]-Crowding_Set[(m+2),1]
          for (i in 2:(ncol(Crowding_Set)-1)){
  	         Crowding_Set[nrow(Crowding_Set),i]<- Crowding_Set[nrow(Crowding_Set),i]+(Crowding_Set[(m+2),(i+1)]-Crowding_Set[(m+2),(i-1)])/Spannweite
          }

        # Setze Distanzmass bzgl. m. Zielfunktion mit kleinstem und groesstem Wert auf Unendlich
          Crowding_Set[nrow(Crowding_Set), which.min(Crowding_Set[(m+2),])]<- Inf
          Crowding_Set[nrow(Crowding_Set), which.max(Crowding_Set[(m+2),])]<- Inf

      Distance_overall<- Distance_overall + Crowding_Set[nrow(Crowding_Set),]
      } # Ende der for-Schleife bzgl. der Berechnung der Crowding-Distanzen der m Zielfunktionen

  # Loesungen mit kleinerem Distanzmass befinden sich in gedraengteren Umgebungen und sind weniger gut!
    Crowding_Set     <- rbind(Crowding_Set, Distance_overall)
    iii              <- order(Crowding_Set[nrow(Crowding_Set),], decreasing = TRUE)
    Crowding_Set     <- Crowding_Set[,iii]
    Crowding_Set_kept<- Crowding_Set[1:4,(1:n_keep)]

  } # Ende der Berechnung des neuen Crowding-Sets, falls das Original aus mehr als 2 Individuen bestand!

  # Neue Population
    Fitness_sorted_Original<- cbind(Fitness_sorted_Original[,which(Fitness_sorted_Original[1,]<cd_rank)], Crowding_Set_kept)

  } # Falls n_toomany > 0, d.h. die Crowding-Distance muss berechnet werden

  # Ausgabe:
    return(Fitness_sorted_Original)
}
