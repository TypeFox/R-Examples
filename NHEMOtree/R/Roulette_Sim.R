Roulette_Sim <-
function(Trees){
  N_Parents<- length(Trees)

  # Rausschreiben der einzelnen Fitnesswerte
    FKR<- c(); Finanzen<- c()
    for (i in 1:N_Parents){
         FKR[i]     <- Trees[[i]]$FKR
         Finanzen[i]<- Trees[[i]]$Finanzen
    }

  # Erstellen eines Fitnesswertes je Individuum
    Temp        <- rbind(1:N_Parents, FKR, Finanzen)
    NDS_Rang    <- nds_rank(Temp[2:3,])
    NDS_Rang_neu<- abs(NDS_Rang-max(NDS_Rang)-1)

  # Berechnung der Auswahlwahrscheinlichkeit 
    Wsk <- NDS_Rang_neu/sum(NDS_Rang_neu) 
    Temp<- rbind(Wsk, Temp)

  # Sortieren der Fitnesswerte
    ii     <- order(Temp[1,])
    Fitness<- Temp[,ii]

  # Selektionsschritt
    E2<- E1<- NULL
    # 1. Elternteil
      j    <- 1
      r    <- runif(n=1, min=0, max=1)
      Summe<- 0
      while(r >= Summe && j<=N_Parents){
            Summe<- Summe+Fitness[1,j]
            E1   <- (Fitness[2,j])
            j    <- j+1
      }

    # 2. Elternteil
      E2 <- E1_temp<- E1
      Anz<- 1
      while(E1_temp==E2){
        j    <- 1
        r    <- runif(n=1, min=0, max=1)
        Summe<- 0
        while(r >= Summe && j<=N_Parents){
              Summe<- Summe+Fitness[1,j]
              E2   <- (Fitness[2,j])
              j    <- j+1; j
        }
      Anz<- Anz+1
      if(Anz>5) E1_temp<- 0  # Zusaetzliches Abbruchkriterium
      }

  # Ausgabe:
    return(c(E1, E2))
}
