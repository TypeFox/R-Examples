Mutation_Punkt_voll <-
function(Tree, N_Varis, uS, oS){
  # Neue Variable und neuer Cutoff-Wert
    Var_neu<- Var_ZHT_Sim(Knoten=1, N_Varis=N_Varis)[1,2]
    CO_neu <- Cutoff_ZHT_Sim(Knoten=1, uS=uS, oS=oS)[1,2]

  # Falls Baum ein Wurzelbaum ist, ersetze Variable in Wurzelknoten
    if (Tree$Knoten==1){
        Tree$Varis[1,2]<- Var_neu
        Tree$CO[1,2]   <- CO_neu
    }

  # Falls Baum kein Wurzelbaum
    if (Tree$Knoten>1){
        Knoten    <- sample(Tree$Varis[,2], 1) # Auswahl des zu mutierenden Knotens
        Knoten_Pos<- which(Tree$Varis[,2]==Knoten)
 
        Tree$Varis[Knoten_Pos,2]<- Var_neu
        Tree$CO[Knoten_Pos,2   ]<- CO_neu
    }

  # Ausgabe:
    return(Tree)
}
