Mutation_overall <-
function(Tree, N_Varis, uS, oS, Max_Knoten, Mutationswsk, Init_Wsk){
  Tree_Mut<- Tree

# Soll mutiert werden?
  temp<- runif(n=1, min=-0.01, max=1.01)

  if (temp<=Mutationswsk){
    # 1) Mutation_Punkt
    # 2) Mutation_Vertauschung
    # 3) Mutation_Permutation 
    # 4) Mutation_Hoist
    # 5) Mutation_Expansion_NEU
    # 6) Delete_Branch_random
    # 7) MBf_max_Knoten

    # Zufaellige Auswahl einer der 7 Mutationsarten
      mb_sample<- sample(1:7, 1)

    # Falls zu mutierender Baum ein Wurzelbaum ist,
      # Punkt-Permutation oder 
      # Vertauschung oder
      # Expansionsmutation oder
      # Subbaum-Mutation

    if (Tree_Mut$Knoten==1) mb_sample<- sample(c(1,2,5,7),1)

    # 1) Punkt-Mutation
      if (mb_sample==1) Tree_Mut<- Mutation_Punkt(Tree=Tree_Mut, N_Varis=N_Varis, uS=uS, oS=oS)

    # 2) Vertauschung
      if (mb_sample==2) Tree_Mut<- Mutation_Vertauschung(Tree=Tree_Mut, N_Varis=N_Varis, uS=uS, oS=oS)

    # 3) Permutation
      if (mb_sample==3) Tree_Mut<- Mutation_Permutation(Tree=Tree_Mut)

    # 4) Hoist-Mutation
      if (mb_sample==4) Tree_Mut<- Mutation_Hoist(Tree=Tree_Mut)

    # 5) Expansionsmutation
      if (mb_sample==5){
          KN          <- Max_Knoten-Tree_Mut$Knoten
          if (KN<0) KN<- 0
          Tree_Mut<- Mutation_Expansion_NEU(Tree=Tree_Mut, Knoten_max_neu=KN, WSK=Init_Wsk, 
                                            N_Varis=N_Varis, uS=uS, oS=oS)
      }  

    # 6) Ast-Loeschung
      if (mb_sample==6) Tree_Mut<- Delete_Branch_random(Tree=Tree_Mut)  

    # 7) Subbaum-Mutation
      if (mb_sample==7){
          Tree_Mut<- MBf_max_Knoten(Tree=Tree_Mut, N_Varis=N_Varis, Knoten_max=Max_Knoten, 
                                    WSK=Init_Wsk, uS=uS, oS=oS)
      }

    } # Ende Mutation

  # Ausgabe:
    return(Tree_Mut)
}
