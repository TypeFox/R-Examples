Best_Tree_new <-
function(data, KOSTEN, Ind, CV){

  # Ermittlung der Fitnesswerte
    FKR              <- sapply(Ind, "[[", "Misclassification")
    Kosten           <- sapply(Ind, "[[", "Costs")
    Finale_Population<- cbind(FKR, Kosten)

  # Baum mit der geringsten FKR
    Finale_Population<- cbind(1:nrow(Finale_Population), Finale_Population)
    ii               <- order(Finale_Population[,2], Finale_Population[,3], decreasing = F)
    Finale_Population<- Finale_Population[ii,]
    Individuum_ID    <- Finale_Population[1,1]
    Variablennamen   <- names(data)[-1]
    proteins_used    <- Chromo_Var(Chromosom=Ind[[Individuum_ID]]$Vector, Gesamtvariablen=Variablennamen)
 
  # Graphik und Ausgabe
    Best_Tree_res<- Ctree_Sim(Data=data, Varis=proteins_used, CV.Laeufe=CV)
    Ausgabe      <- list()
    Ausgabe[[1]] <- Best_Tree_res[[6]]                                         # rpart objects (pruned tree)
    Ausgabe[[2]] <- Best_Tree_res[[4]]                                         # Misclassification rate
    Ausgabe[[3]] <- Cost_calculation(Varis=Best_Tree_res[[3]], KOSTEN=KOSTEN)  # Costs of pruned tree
    return(Ausgabe)
}
