PF_1perm <-
function(Tree, VIM_alt, Daten_Perm, Daten, CV_Laeufe){
  # Vorbereitung
    VIM          <- VIM_alt
    Spalten_in_DB<- Tree$Varis[,2]
    Spalten_in_DB<- NoDupKey(Spalten_in_DB)

  # Fitness des Originalbaumes
    P_N         <- Spalten_in_DB-1
    P_N         <- NoDupKey(P_N)
    FKR_Original<- 100*FKR_CV(Tree=Tree, Daten=Daten, CV_Laeufe=CV_Laeufe)

  # Fuer alle Variablen Berechnung der Permutierten Fehlerfreiheit
    for (i in 1:length(Spalten_in_DB)){
       # Permutation
         Var_Perm                   <- Daten
         Var_Perm[,Spalten_in_DB[i]]<- Daten_Perm[,Spalten_in_DB[i]]

       # FKR- und VIM-Berechnung
         FKR_Perm<- 100*FKR_CV(Tree=Tree, Daten=Var_Perm, CV_Laeufe=CV_Laeufe)
         FKR_Diss<- FKR_Perm-FKR_Original
         if (FKR_Diss<0) FKR_Diss=0
       
       VIM[Spalten_in_DB[i]-1]<- VIM[Spalten_in_DB[i]-1]+FKR_Diss
    }

  # Ausgabe:
    return(VIM)
}
