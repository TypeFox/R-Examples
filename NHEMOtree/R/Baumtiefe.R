Baumtiefe <-
function(Tree){
  ZHT  <- Tree$Zusammenhangstabelle 
  Tiefe<- rep(0, nrow(ZHT))
  ZHT2 <- cbind(ZHT, Tiefe)

  # Initialisierung:
    i<- 1; temp_i<- 0
    ZHT2[1:2, 4]<- 1

  # Fuer welche Knoten hat der Vaterknoten die Tiefe i: 
    while (temp_i[1]!=-99){
           temp_i<- ZHT2[which(ZHT2[,4]==i), 3]
           for (j in 1:length(temp_i)) ZHT2[which(ZHT2[,1]==temp_i[j]), 4]<- i+1
           i<- i+1
           temp_i<- sort(temp_i, decreasing = TRUE)
    }

  # Ausgabe:
    Tree$Zusammenhangstabelle<- ZHT2
    return(Tree)
}
