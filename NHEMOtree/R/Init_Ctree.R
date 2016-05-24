Init_Ctree <-
function(Daten, KOSTEN, popsize, 
                      in_percent, Chromo_length, 
                      CV.Lauf, proteins){
  Individuen<- list()
  Max_Costs <- sum(KOSTEN[,2])
  Max_Anz   <- nrow(KOSTEN)

  for (i in 1:popsize){
  Tupel_Init<- list()

  # Stochastische Initialisierung
    Chromo_Parent<- Initialisierung_stoch(pop_size=1, n=Chromo_length, in_percent=in_percent)[[1]]
    while (sum(Chromo_Parent[[1]])==0){  # Keine Individuen ohne Variablen
               Chromo_Parent<- Initialisierung_stoch(pop_size=1, n=Chromo_length, in_percent=in_percent)[[1]]
    }

  # Berechne Fitness fuer ein Individuum
    Tupel_Init<- Fitness_Calc_Ctree_Sim(Daten=Daten, CV.Lauf=CV.Lauf, Chromosom=Chromo_Parent[[1]], proteins=proteins)

  # Individuum
    Individuum<- list()

  # Kosten im Pruned-Baum    

    # Umschreiben der Pruned-Variablen in 0/1-Vektoren     
      VectorOut_P <- Chromo_bit_long(Varis=Tupel_Init$P_Tree, Gesamtvariablen=proteins)
      ProteinOut_P<- which(VectorOut_P[[1]]==1)

    # Falls Individuum keine Variablen enthaelt, maximale FKR und Kosten
      FKR_P         <- 100
      CostsOut_P_sum<- Max_Costs
      CostsOut_P_Anz<- Max_Anz

      if (sum(VectorOut_P[[1]])>0){
          FKR_P         <- Fitness_Calc_Ctree_Sim(Daten=Daten, CV.Lauf=CV.Lauf, Chromosom=VectorOut_P[[1]], proteins=proteins)$Misclassification      
          CostsOut_P    <- KOSTEN[ProteinOut_P, 2]
          CostsOut_P_sum<- sum(CostsOut_P)
          CostsOut_P_Anz<- length(CostsOut_P)
      }

    # Individuumerstellung          
      Individuum$Vector           <- VectorOut_P[[1]]
      Individuum$Misclassification<- FKR_P
      Individuum$Costs            <- CostsOut_P_sum
      Individuum$Anz              <- CostsOut_P_Anz

    Individuen[[i]]<- Individuum
  } # Schleife zu Ende fuer Initial-Population

  # Ausgabe
    return(Individuen)
}
