Fitness_Calc_Ctree_Sim <-
function(Daten, CV.Lauf, Chromosom, proteins){
  Tupel_Init<- list()
 
  # Ausgabe der Variablennamen pro enthaltendem Gen:
    proteins_used<- Chromo_Var(Chromosom=Chromosom, Gesamtvariablen=proteins)
   
  # Keine Variablen enthalten
    if (length(proteins_used)==0){
        Tupel_Init              <- list()
        Tupel_Init$Modell       <- paste(names(Daten)[1], "~")
        Tupel_Init$Used_Proteins<- NULL
        Tupel_Init$O_Tree       <- NULL
        Tupel_Init$P_Tree       <- NULL
        Tupel_Init$Misclassification<- 100  #table(Daten[,1])
        return(Tupel_Init)
    }

  # Berechnung Klassifikation:
    Ctree            <- Ctree_Sim(Data=Daten, Varis=proteins_used, CV.Laeufe=CV.Lauf)
    Modell           <- Ctree$Modell
    Original_tree    <- Ctree$Originalbaum
    Pruned_tree      <- Ctree$PrunedBaum
    Misclassification<- min(Ctree$FKR)

  # Schreiben der Ergebnisse hinter der der Eltern:
    Tupel_Init<- list(Modell=Modell, Used_Proteins=proteins_used, O_Tree=Original_tree, P_Tree=Pruned_tree, Misclassification=Misclassification) 

  # Ausgabe
    return(Tupel_Init)
}
