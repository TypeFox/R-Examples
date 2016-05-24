Var_Rang <-
function(Tree, Knoten, VIM, N_Varis){
  Knoten<- NoDupKey(Knoten)
  Knoten_b<- c()

  for (p in 1:length(Knoten)) Knoten_b<- c(Knoten_b, which(Tree$Varis[,1]==Knoten[p]))
  Var<- Tree$Varis[Knoten_b,2]-1

  # Wahrscheinlichkeitsberechnung auf den Raengen
    if (length(VIM)==1) VIM<- rep(1,N_Varis)
    Wichtigkeit<- rank(VIM[Var], ties.method = "average")

    N  <- length(Wichtigkeit)
    Wsk<- c()
    for (i in 1:N) Wsk<- c(Wsk, rep(Knoten_b[i], Wichtigkeit[i]))

  # Ausgabe der Position des Rekombinationsknotens:
    if (length(Wsk)==1) return(Wsk)
    if (length(Wsk)>1) return(sample(Wsk, 1))
}
