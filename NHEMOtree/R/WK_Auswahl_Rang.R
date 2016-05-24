WK_Auswahl_Rang <-
function(Tree, VIM, N_Varis){
  # Abbruch, falls Baum ein Wurzelbaum ist.
    if (Tree$Knoten==1) stop(paste("'WK_Auswahl_Rang' does not allow a tree with just a root node!"))

  # Variablen-Extraktion ohne Variable des Wurzelknotens
    Var<- Tree$Varis[,2][-1]
    Var<- Var-1

  # Wahrscheinlichkeitsberechnung auf den Raengen
    if (length(VIM)==1) VIM<- rep(1,N_Varis)
    Wichtigkeit<- rank(VIM[Var], ties.method = "average")
    N          <- length(Wichtigkeit)
    Wsk        <- c()
    for (i in 1:N) Wsk<- c(Wsk, rep(Var[i], Wichtigkeit[i]))

  # Ausgabe des Rekombinationsknotens:
    if (length(Wsk)==1) Variable<- Wsk+1 # Variablenbezeichnung im Baum

    if (length(Wsk)>1){
        Variable<- sample(Wsk, 1)
        Variable<- Variable+1            # Variablenbezeichnung im Baum
    }

  # Falls Variablen mehrfach im Baum vorkommt:
    Knoten<- which(Tree$Varis[,2]==Variable)
    if (Knoten[1]==1)       Knoten<- Knoten[-1]
    if (length(Knoten) > 1) Knoten<- sample(Knoten,1)

  # Ausgabe:
    return(Knoten)
}
