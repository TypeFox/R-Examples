Init_RHaH_Sim <-
function(Daten, N_Varis, Proteinkosten, 
                         Popsize, Max_Knoten=Max_Knoten, CV.Lauf, 
                         Init_Wsk=Init_Wsk, uS=uS, oS=oS){

    Baum<- Init_RHaH(N=Popsize, Knoten_max=Max_Knoten, Wsk=Init_Wsk, N_Varis=N_Varis, uS=uS, oS=oS)

    Rate<- c(); Finanzen<- c(); Anzahl<- c(); Knoten<- c()
    Quantile_Knoten<- NULL
    Quantile_FKR   <- NULL

    for (i in 1:Popsize){
      # Misclassification
        Rate[i]<- 100*FKR_CV(Tree=Baum[[i]], Daten=Daten, CV_Laeufe=CV.Lauf)

      # Costs
        P1<- Baum[[i]]$Varis[,2]-1 # Protein-IDs (ohne Datenbank)
        P2<- NoDupKey(P1)
      
        Knoten[i]  <- Baum[[i]]$Knoten
        Finanzen[i]<- getCosts_Tree_sum_Sim(Kostenmatrix=Proteinkosten, Protein_IDs=P2)
        Anzahl[i]  <- length(P2)

      # Tree
        Baum[[i]]$FKR     <- Rate[i]   
        Baum[[i]]$Finanzen<- Finanzen[i]
        Baum[[i]]$Anzahl  <- Anzahl[i]
    } 

  # Output
    Ausgabe         <- list()
    Ausgabe$Knoten  <- Knoten
    Ausgabe$Rate    <- Rate
    Ausgabe$Finanzen<- Finanzen
    Ausgabe$Anzahl  <- Anzahl
    Ausgabe$Baum    <- Baum
    return(Ausgabe)
}
