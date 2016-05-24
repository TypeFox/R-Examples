X_Brood_VIM2 <-
function(Daten, E1, E2, N_Varis, uS, oS,
                        Max_Knoten, Brutgroesse, Kostenmatrix, CV_Laeufe, VIM){

  if (E1$Knoten==1 && E2$Knoten==1) return(list(E1,E2))
  
  # Kindererstellung:
    Temp       <- list()
    Kinder     <- list()
    Kinder_Best<- list()

    # Falls Eltern keine Wurzelbaeume:
      if (E1$Knoten!=1 && E2$Knoten!=1){
          for (i in 1:Brutgroesse){
               Temp           <- X_SGP_VIM2(E1=E1, E2=E2, 
                                            Max_Knoten=Max_Knoten, 
                                            Daten=Daten, CV_Laeufe=CV_Laeufe, VIM=VIM,
                                            N_Varis=N_Varis, uS=uS, oS=oS)

               Kinder[[2*i-1]]<- Temp[[1]]
               Kinder[[2*i]]  <- Temp[[2]]    
      }}

    # Falls ein Elternteil ein Wurzelbaum:
      if (E1$Knoten==1 || E2$Knoten==1){
          for (i in 1:Brutgroesse){   
               Temp           <- X_SGP_WK1_VIM(E1=E1, E2=E2, VIM=VIM, N_Varis=N_Varis)
               Kinder[[2*i-1]]<- Temp[[1]]
               Kinder[[2*i]]  <- Temp[[2]]  
      }}

  # Fitnesserstellung:
    Fitness_temp<- Fitness_Nachkommen_Sim(Trees=Kinder, Daten=Daten, Kostenmatrix=Kostenmatrix, 
                                          CV_Laeufe=CV_Laeufe)
    Fitness_temp<- rbind(1:(2*Brutgroesse), Fitness_temp)

  # Vergleich der Fitnesswerte
    Fitness       <- Fitness_temp[-4,]
    Fitness       <- rbind(nds_rank(Fitness[2:3,]), Fitness)
    ii            <- order(Fitness[1,], Fitness[2,], Fitness[3,], Fitness[4,])
    Fitness_sorted<- Fitness[,ii]

  # Rang an Position 2
    cd_rank      <- Fitness_sorted[1,2]
    which_cd_rank<- which(Fitness_sorted[1,]<=cd_rank)
    n_toomany    <- length(which_cd_rank)-2

  # Falls Crowding Distance Sorting NICHT notwendig
    if (n_toomany==0){
        Fitness_sorted<- Fitness_sorted[,which(Fitness_sorted[1,]<=cd_rank)]
    }

  # Falls Crowding Distance Sorting notwendig
    if (n_toomany > 0){
      # Alle Individuen mit Rang <= 'Temp'-Rang
        FS_after_NDS<- Fitness_sorted[,which_cd_rank]

      # Berechnung der Crowding-Distance fuer alle Individuen mit Rang 'cd_rank'
        Crowding_Set<- FS_after_NDS[,which(FS_after_NDS[1,]==cd_rank)]
        n_keep      <- ncol(Crowding_Set)-n_toomany

      # Falls Crowding-Set nur aus zwei Individuen besteht --> Zufall
        if (ncol(Crowding_Set)==2){
            temp<- runif(1,0,1)
          if (temp <  0.5) Crowding_Set_kept<- Crowding_Set[,1]
          if (temp >= 0.5) Crowding_Set_kept<- Crowding_Set[,2]
        }

        if (ncol(Crowding_Set)>2){
        # Berechne Distanzmasse fuer die m Zielfunktionen
          n_objectives    <- nrow(Crowding_Set)-2     # Anzahl der Zielfunktionen
          Distance_overall<- rep(0, ncol(Crowding_Set))

          for (m in 1:n_objectives){
          # Sortieren
            ii          <- order(Crowding_Set[(m+2),])
            Crowding_Set<- Crowding_Set[,ii]
            Dim         <- dim(Crowding_Set)
            Distance    <- rep(0, Dim[2])
            Crowding_Set<- rbind(Crowding_Set, Distance)

          # Berechne Distanzmass bzgl. m. Zielfunktion fuer alle Werte, die keine Extremwerte sind.
            Spannweite<- Crowding_Set[(m+2),Dim[2]]-Crowding_Set[(m+2),1]
            for (i in 2:(Dim[2]-1)){
                   Crowding_Set[(Dim[1]+1),i]<- Crowding_Set[(Dim[1]+1),i]+(Crowding_Set[(m+2),(i+1)]-Crowding_Set[(m+2),(i-1)])/Spannweite
            }
            
          # Setze Distanzmass bzgl. m. Zielfunktion mit kleinstem Wert auf Inf
            Crowding_Set[(Dim[1]+1), which.min(Crowding_Set[(m+2),])]<- Inf

          # Setze Distanzmass bzgl. m. Zielfunktion mit groesstem Wert auf Inf
            Crowding_Set[(Dim[1]+1), which.max(Crowding_Set[(m+2),])]<- Inf

          Distance_overall<- Distance_overall + Crowding_Set[(Dim[1]+1),]
          } # Ende der for-Schleife bzgl. der Berechnung der Crowding-Distanzen 
            # der m Zielfunktionen

        # Loesungen mit kleinerem Distanzmass befinden sich in gedraengteren 
        # Umgebungen und sind weniger gut!
          Crowding_Set     <- rbind(Crowding_Set, Distance_overall)
          iii              <- order(Crowding_Set[nrow(Crowding_Set),], decreasing = TRUE)
          Crowding_Set     <- Crowding_Set[,iii]
          Crowding_Set_kept<- Crowding_Set[1:(2+n_objectives),(1:n_keep)]
        } # Ende der Berechnung des neuen Crowding-Sets, falls das Original aus mehr als 2 Individuen bestand!
  
      # Neue Population
        Fitness_sorted<- cbind(Fitness_sorted[,which(Fitness_sorted[1,]<cd_rank)], Crowding_Set_kept)

      } # Falls n_toomany > 0, d.h. die Crowding-Distance muss berechnet werden

  # ENDE - Vergleich der Fitnesswerte

  # Die besten Kinder:
    Kinder_Best[[1]]<- Kinder[[Fitness_sorted[2,1]]]
    Kinder_Best[[2]]<- Kinder[[Fitness_sorted[2,2]]]

  # Einfuegen der Fitnesswerte
    Kinder_Best[[1]]$FKR     <- Fitness_temp[2,Fitness_sorted[2,1]]
    Kinder_Best[[2]]$FKR     <- Fitness_temp[2,Fitness_sorted[2,2]]

    Kinder_Best[[1]]$Finanzen<- Fitness_temp[3,Fitness_sorted[2,1]]
    Kinder_Best[[2]]$Finanzen<- Fitness_temp[3,Fitness_sorted[2,2]]

    Kinder_Best[[1]]$Anzahl  <- Fitness_temp[4,Fitness_sorted[2,1]]
    Kinder_Best[[2]]$Anzahl  <- Fitness_temp[4,Fitness_sorted[2,2]]

  # Ausgabe der besten Kinder samt Fitnesswerte
    return(Kinder_Best)
}
