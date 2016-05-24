Mutation_Permutation <-
function(Tree){
  # Falls Tree Wurzelbaum, keine Permutation moeglich:
    if (Tree$Knoten==1) return(Tree)

  # Auswahl eines Permutationsknotens mit weniger als zwei Blaetter als Tochterknoten
    ZHT  <- Tree$Zusammenhangstabelle
    temp <- ZHT[which(ZHT[,3]!=-99), 1]

    # Falls der Baum nur einen Tochterknoten hat, vertausche Tochterknoten mit Blatt auf der anderen Seite des Wurzelknoten:
      if (length(temp)==1){
          Zeile<- which(ZHT[,3]!=-99)
          ZHT_P_sorted<- ZHT
          if (ZHT[Zeile, 2]==0){  # Tochterknoten links vom Wurzelknoten
              ZHT_P_sorted[1,3]<- -99
              ZHT_P_sorted[2,3]<- ZHT[Zeile,3]
          }
          if (ZHT[Zeile, 2]==1){  # Tochterknoten rechts vom Wurzelknoten
              ZHT_P_sorted[1,3]<- ZHT[Zeile,3]
              ZHT_P_sorted[2,3]<- -99
          }
       }

    # Falls der Baum mehr als einen Tochterknoten hat, vertausche nur die Tochterknoten:
      if (length(temp)>1){ 
          temp2<- temp[1]
  
          for (i in 2:length(temp)){
               if (temp[i]!=temp[i-1]) temp2<- c(temp2, temp[i])
          }

        # Auswahl des Knotens, dessen Toechterknoten permutiert werden sollen    
          P_temp<- which(ZHT[,1]==sample(temp2, 1))

          # Vertauschen der Positionen der Unterbaeume des Permutationsknotens
            ZHT_P             <- ZHT
            ZHT_P[P_temp[1],2]<- 1
            ZHT_P[P_temp[2],2]<- 0
 
          # Sortieren der ZHT, so dass stets zuerst Position 0 und dann Position 1
            ii    <- order(ZHT_P[,1], ZHT_P[,2], ZHT_P[,3])
            ZHT_P_sorted<- ZHT_P[ii,] 
      }

  # Ausgabe:
    Tree$Zusammenhangstabelle<- ZHT_P_sorted
    return(Tree)
}
