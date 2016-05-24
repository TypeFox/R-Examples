Delete_EL <-
function(Tree, Daten, N_Varis, uS=uS, oS=oS){
  # Bestimmung der leeren Blaetter
    Temp_LB<- Empty_Laeves(Tree=Tree, Daten=Daten)

  # Loeschen
    Tree_new      <- Temp_LB[[1]]
    Leere_Blaetter<- Temp_LB[[2]]
    while (Leere_Blaetter[1] > 1){
        Tree_new      <- Delete_Knot_fixed(Tree=Tree_new, Knoten=Leere_Blaetter[1], N_Varis=N_Varis, uS=uS, oS=oS)
        Temp_LB       <- Empty_Laeves(Tree=Tree_new, Daten=Daten)
        Tree_new      <- Temp_LB[[1]]
        Leere_Blaetter<- Temp_LB[[2]]
    }
    if (Leere_Blaetter[1] == 1) Tree_new<- Delete_Knot_fixed(Tree=Tree_new, Knoten=Leere_Blaetter[1], 
                                                             N_Varis=N_Varis, uS=uS, oS=oS) 

  # Ausgabe:
    return(Tree_new)
}
