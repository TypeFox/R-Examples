Mutation_Punkt_Vari <-
function(Tree, N_Varis){
  # Neue Variable
    Var_neu<- Var_ZHT_Sim(Knoten=1, N_Varis=N_Varis)[1,2]

  # Falls Baum ein Wurzelbaum ist, ersetze Variable in Wurzelknoten
    if (Tree$Knoten==1) Tree$Varis[1,2]<- Var_neu

  # Falls Baum kein Wurzelbaum
    if (Tree$Knoten>1){
        Knoten    <- sample(Tree$Varis[,2], 1) 
        Knoten_Pos<- which(Tree$Varis[,2]==Knoten)
 
        Tree$Varis[Knoten_Pos,2]<- Var_neu
    }

  # Ausgabe:
    return(Tree)
}
