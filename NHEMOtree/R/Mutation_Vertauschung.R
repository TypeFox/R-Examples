Mutation_Vertauschung <-
function(Tree, N_Varis, uS, oS){
  # Falls Baum ein Wurzelbaum ist, erstelle einen neuen Wurzelbaum
    if (Tree$Knoten==1) Tree<- Baum_zufaellig_komplett_Sim(Knoten_max=1, N_Varis=N_Varis, uS=uS, oS=oS)

  # Falls Baum kein Wurzelbaum
    if (Tree$Knoten>1){
        Knoten<- sample(1:Tree$Knoten, 2)
        
        Vari1 <- Tree$Varis[Knoten[1],2]
        Vari2 <- Tree$Varis[Knoten[2],2]
        CO1   <- Tree$CO[Knoten[1],2]
        CO2   <- Tree$CO[Knoten[2],2]

        # Verauschen
          Tree$Varis[Knoten[1],2]<- Vari2
          Tree$Varis[Knoten[2],2]<- Vari1
          Tree$CO[Knoten[1],2]   <- CO2
          Tree$CO[Knoten[2],2]   <- CO1
    }

  # Ausgabe:
    return(Tree)
}
