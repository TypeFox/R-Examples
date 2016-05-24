Mutation_Punkt <-
function(Tree, N_Varis, uS, oS){
  temp<- sample(1:2,1)

  # Mutation nur der Variablen
    if (temp==1) Ausgabe<- Mutation_Punkt_Vari(Tree=Tree, N_Varis=N_Varis)

  # Mutation der Variablen und des Cutoffs
    if (temp==2) Ausgabe<- Mutation_Punkt_voll(Tree=Tree, N_Varis=N_Varis, uS=uS, oS=oS)

  # Ausgabe:
    return(Ausgabe)
}
