TS_Sim <-
function(Tree, K=4){
  N_Parents<- length(Tree)
  if (N_Parents < K) stop(paste("Tournament size too big!"))

  Eltern1<- sample(x=1:N_Parents, size=K, replace = FALSE)
  Eltern2<- sample(x=1:N_Parents, size=K, replace = FALSE)

  ################
  # Elternteil 1 #
  ################     
    FKR1<- c(); Finanzen1<- c()
    for (i in 1:K) FKR1[i]     <- Tree[[Eltern1[i]]]$FKR 
    for (i in 1:K) Finanzen1[i]<- Tree[[Eltern1[i]]]$Finanzen

    Temp1     <- rbind(Eltern1, FKR1, Finanzen1)
    Vergleich1<- nds_rank(Temp1[2:3,])
    Min1      <- which(Vergleich1==min(Vergleich1))

    if (length(Min1)==1) E1<- Eltern1[Min1]  # Ein Elternteil mit kleinstem nds_rank

    if (length(Min1)==2){ # Zwei Elternteile mit kleinstem nds_rank --> Zufall entscheidet
        temp<- runif(1,0,1)
        if (temp <  0.5) E1<- Eltern1[Min1[1]]
        if (temp >= 0.5) E1<- Eltern1[Min1[2]]
    }
 
    if (length(Min1)>2){  # Crowding-Distanz kommt zum Tragen
        Crowding_Set1   <- Temp1[,Min1]
        Distance_overall<- matrix(rep(0, ncol(Crowding_Set1)), ncol=1)

        for (m in 1:2){  # 2 Zielfunktionen werden betrachtet!
        # Sortieren
          ii              <- order(Crowding_Set1[(m+1),])
          Crowding_Set1   <- Crowding_Set1[,ii]
          Distance_overall<- t(as.matrix(Distance_overall))[,ii]
          Distance        <- rep(0, ncol(Crowding_Set1))
          Crowding_Set1   <- rbind(Crowding_Set1, Distance)

        # Berechne Distanzmass bzgl. m. Zielfunktion fuer alle Werte, die keine Extremwerte sind.
          Spannweite<- Crowding_Set1[(m+1), ncol(Crowding_Set1)]-Crowding_Set1[(m+1),1]
          for (i in 2:(ncol(Crowding_Set1)-1)){
               Crowding_Set1[nrow(Crowding_Set1),i]<- Crowding_Set1[nrow(Crowding_Set1),i]+(Crowding_Set1[(m+1),(i+1)]-Crowding_Set1[(m+1),(i-1)])/Spannweite
          }

        # Setze Distanzmass bzgl. m. Zielfunktion mit kleinstem Wert auf Unendlich
          Crowding_Set1[nrow(Crowding_Set1), which.min(Crowding_Set1[(m+1),])]<- Inf

        # Setze Distanzmass bzgl. m. Zielfunktion mit groesstem Wert auf Unendlich
          Crowding_Set1[nrow(Crowding_Set1), which.max(Crowding_Set1[(m+1),])]<- Inf

          Distance_overall<- Distance_overall + Crowding_Set1[nrow(Crowding_Set1),]
      } # Ende der for-Schleife bzgl. der Berechnung der Crowding-Distanzen der m Zielfunktionen

      # Loesungen mit kleinerem Distanzmass befinden sich in gedraengteren Umgebungen und sind weniger gut!
      # --> Auswahl des Individuums mit der groessten Crowding-Distanz
        Crowding_Set1<- rbind(Crowding_Set1, Distance_overall)
        E1           <- Crowding_Set1[1, sample(x=which(Distance_overall==Inf), size=1)]
      } 

  ################
  # Elternteil 2 #
  ################
    E2 <- E1_temp<- E1
    Anz<- 1
    while(E1_temp==E2){
    FKR2<- c(); Finanzen2<- c()
    for (i in 1:K) FKR2[i]     <- Tree[[Eltern2[i]]]$FKR 
    for (i in 1:K) Finanzen2[i]<- Tree[[Eltern2[i]]]$Finanzen

    Temp2     <- rbind(Eltern2, FKR2, Finanzen2)
    Vergleich2<- nds_rank(Temp2[2:3,])
    Min2      <- which(Vergleich2==min(Vergleich2))
    
    if (length(Min2)==1) E2<- Eltern2[Min2]  # Ein Elternteil mit kleinstem nds_rank
 
    if (length(Min2)==2){ # Zwei Elternteile mit kleinstem nds_rank --> Zufall entscheidet
        temp<- runif(1,0,1)
        if (temp <  0.5) E2<- Eltern2[Min2[1]]
        if (temp >= 0.5) E2<- Eltern2[Min2[2]]
    }

    if (length(Min2)>2){  # Crowding-Distanz kommt zum Tragen
        Crowding_Set2   <- Temp2[,Min2]
        Distance_overall<- matrix(rep(0, ncol(Crowding_Set2)), ncol=1)

        for (m in 1:2){  # 2 Zielfunktionen werden betrachtet!
        # Sortieren
          ii              <- order(Crowding_Set2[(m+1),])
          Crowding_Set2   <- Crowding_Set2[,ii]
          Distance_overall<- t(as.matrix(Distance_overall))[,ii]
          Distance        <- rep(0, ncol(Crowding_Set2))
          Crowding_Set2   <- rbind(Crowding_Set2, Distance)

        # Setze Distanzmass bzgl. m. Zielfunktion mit kleinstem Wert auf Unendlich
          Crowding_Set2[nrow(Crowding_Set2), which.min(Crowding_Set2[(m+1),])]<- Inf

        # Setze Distanzmass bzgl. m. Zielfunktion mit groesstem Wert auf Unendlich
          Crowding_Set2[nrow(Crowding_Set2), which.max(Crowding_Set2[(m+1),])]<- Inf

        # Berechne Distanzmass bzgl. m. Zielfunktion fuer alle Werte, die keine Extremwerte sind.
          Spannweite<- Crowding_Set2[(m+1), ncol(Crowding_Set2)]-Crowding_Set2[(m+1),1]
          for (i in 2:(ncol(Crowding_Set2)-1)){
               Crowding_Set2[nrow(Crowding_Set2),i]<- Crowding_Set2[nrow(Crowding_Set2),i]+(Crowding_Set2[(m+1),(i+1)]-Crowding_Set2[(m+1),(i-1)])/Spannweite
          }

          Distance_overall<- Distance_overall + Crowding_Set2[nrow(Crowding_Set2),]
        } # Ende 'for (m in 1:2)'

      # Loesungen mit kleinerem Distanzmass befinden sich in gedraengteren Umgebungen und sind weniger gut!
      # --> Auswahl des Individuums mit der groessten Crowding-Distanz
        Crowding_Set2<- rbind(Crowding_Set2, Distance_overall)
        E2           <- Crowding_Set2[1, sample(x=which(Distance_overall==Inf), size=1)]
      } 
      Anz<- Anz+1; # print(Anz)
      if(Anz>5) E1_temp<- 0  # Zusaetzliches Abbruchkriterium
    }

  # Ausgabe der Eltern
    return(c(E1, E2))
}
