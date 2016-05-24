MBf_max_Knoten <-
function(Tree, N_Varis, Knoten_max, WSK, Wsk_Abnahme=0.1, uS, oS){

  # Einmaliges Ersetzen eines Astes   
    Baum_Mutiert<- Mutation_Branch_fully_Sim(Tree=Tree, N_Varis=N_Varis, Knoten_max=Knoten_max, 
                                             WSK=WSK, Wsk_Abnahme=Wsk_Abnahme, uS=uS, oS=oS)

  # Falls neuer Baum zu gross ist, ersetze zu lange Aeste bis zur maximalen Groesse
    while (Baum_Mutiert$Knoten>Knoten_max){
      Baum_Mutiert<- Mutation_Branch_fully_Sim(Tree=Baum_Mutiert, N_Varis=N_Varis, Knoten_max=Knoten_max, 
                                               WSK=WSK, Wsk_Abnahme=Wsk_Abnahme, uS=uS, oS=oS)
    }

  # Ausgabe:
    return(Baum_Mutiert)
}
