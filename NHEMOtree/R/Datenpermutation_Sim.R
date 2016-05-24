Datenpermutation_Sim <-
function(Daten, Vari_N){
  Daten_Perm<- Daten

  # Permutation 
    for (i in 1:Vari_N) Daten_Perm[,i+1]<- sample(Daten_Perm[,i+1])

  # Ausgabe:
    return(Daten_Perm)
}
