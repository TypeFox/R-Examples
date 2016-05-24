Init_RHaH_opt_Sim_Gini <-
function(Daten, N_Varis, Proteinkosten, 
                                  Popsize, Max_Knoten, CV.Lauf, 
                                  Init_Wsk, uS, oS, AnzCO){
  # Baumerstellung
    Baum_temp<- Init_RHaH(N=Popsize, Knoten_max=Max_Knoten, Wsk=Init_Wsk, N_Varis=N_Varis, uS=uS, oS=oS)
    Baum<- list(); Rate<- c(); Finanzen<- c(); Anzahl<- c(); Knoten<- c()
    Quantile_Knoten<- NULL
    Quantile_FKR   <- NULL
    
    for (p in 1:Popsize){
      # Bestimmung der optimalen Cutoffs
        Baum[[p]]<- Cutoffwahl_gini_Sim(Tree=Baum_temp[[p]], Daten=Daten, CV.Lauf=CV.Lauf, AnzCO=AnzCO)

      # Fitnessberechnungen
        Rate[p]<- 100*FKR_CV(Tree=Baum[[p]], Daten=Daten, CV_Laeufe=CV.Lauf)
      
      # Kostenberechnungen 
        P1         <- Baum[[p]]$Varis[,2]-1 # Protein-IDs (ohne Datenbank)
        P2         <- NoDupKey(P1)
        Knoten[p]  <- Baum[[p]]$Knoten
        Finanzen[p]<- getCosts_Tree_sum_Sim(Kostenmatrix=Proteinkosten, Protein_IDs=P2)
        Anzahl[p]  <- length(P2)
      
      # Schreiben der Fitnesswerte in 'Baum'
        Baum[[p]]$FKR     <- Rate[p]   
        Baum[[p]]$Finanzen<- Finanzen[p]
        Baum[[p]]$Anzahl  <- Anzahl[p]
    } 
    
    # Vorbereitung der Ausgabe
      Ausgabe         <- list()
      Ausgabe$Knoten  <- Knoten
      Ausgabe$Rate    <- Rate
      Ausgabe$Finanzen<- Finanzen
      Ausgabe$Anzahl  <- Anzahl
      Ausgabe$Baum    <- Baum
    
    # Ausgabe: 
      return(Ausgabe)
}
