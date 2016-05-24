Empty_Laeves <-
function(Tree, Daten){
  Tree     <- Werte_Reduktion_ZHT2(Tree)  # Korrekte Nummerierung der ZHT-Spalteneintraege
  ZHT      <- Tree$Zusammenhangstabelle
  N_Knoten <- Tree$Knoten
  Variablen<- Tree$Varis
  Cutoffs  <- Tree$CO
  BeobNeu  <- list()
  FK       <- c()    # Fehlklassifikation pro Analysepunkt (absolute Zahlen)

  # Falls Baum Wurzelbaum ist:
    if (nrow(ZHT)==2){
        # Welche Datenpunkte erfuellen Bedingung
          BeobNeu[[1]]<- which(Daten[, Variablen[1,2]] <= Cutoffs[1,2]) 
          BeobNeu[[2]]<- which(Daten[, Variablen[1,2]] >  Cutoffs[1,2]) 
    }

  # Falls Baum kein Wurzelbaum ist
    if (nrow(ZHT)>2){
        BeobNeu[[1]]<- which(Daten[, Variablen[1,2]] <= Cutoffs[1,2]) 
        BeobNeu[[2]]<- which(Daten[, Variablen[1,2]] >  Cutoffs[1,2]) 
   
        for(i in 2:N_Knoten){
            r<- ZHT[which(ZHT[,3]==i),]      # Zeilen mit 3.Spalte = i 
            k<- 2*r[1]-(1-r[2])
            BeobNeu[[2*i-1]]<- intersect(which(Daten[, Variablen[i,2]] <= Cutoffs[i,2]), BeobNeu[[k]]) 
            BeobNeu[[2*i]]  <- intersect(which(Daten[, Variablen[i,2]] >  Cutoffs[i,2]), BeobNeu[[k]])
        }
     }

  # Leere Blaetter
    Leaves<- c()
    for (m in 1:length(BeobNeu)) Leaves[m]<- length(BeobNeu[[m]])

    Leaves_empty<- ZHT[which(Leaves==0),1]

    # Sortieren und Entfernen aller Doppelten
      temp1<- sort(Leaves_empty)
      temp2<- temp1[1]
      if (length(temp1)==0) temp2<- 0
      if (length(temp1)>1){
          for (l in 1:(length(temp1)-1)){
               if (temp1[l]!=temp1[l+1]) temp2<- c(temp2, temp1[l+1])
          }  
      }

      Leaves_empty<- sort(temp2, decreasing=TRUE)

  # Ausgabe:
    return(list(Tree=Tree, Leaves_empty=Leaves_empty))
}
