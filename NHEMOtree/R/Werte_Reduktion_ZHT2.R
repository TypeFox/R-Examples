Werte_Reduktion_ZHT2 <-
function(Tree){
  ZHT<- Tree$Zusammenhangstabelle

  # Sortieren der Vaterknoten
    ii     <- order(ZHT[,1], ZHT[,2], ZHT[,3])
    ZHT_neu<- ZHT[ii,]
    ZHT    <- ZHT_neu

  # Neuer Wurzelknoten mit Wert 1:
    Abzug<- ZHT[1,1]-1
    ZHT[,1]<- ZHT[,1]-Abzug 
    for (l in 1:nrow(ZHT)){
         if (ZHT[l,3]!=-99) ZHT[l,3]<-ZHT[l,3]-Abzug
    }

  # Eintraege in der 1. ZHT-Spalte
    ZHT_1 <- ZHT[2*(1:(nrow(ZHT)/2)), 1]
    ZHT_1L<- length(ZHT_1)

  # Reduzierung der Werte in der 1. und 3. ZHT-Spalte
    i<- 1
    while (i < ZHT_1L){
      if ((ZHT_1[i]+1)!=ZHT_1[i+1]){
           Diff<- ZHT_1[i+1]-(ZHT_1[i]+1)

         # Abziehen der Differenz von allen Eintraegen = ZHT_1[i+1] in ZHT
           ZHT[which(ZHT[,1]==ZHT_1[i+1]),1]<- ZHT_1[i+1] - Diff
           ZHT[which(ZHT[,3]==ZHT_1[i+1]),3]<- ZHT_1[i+1] - Diff

         # Erneuern von ZHT_1
           ZHT_1<- ZHT[2*(1:(nrow(ZHT)/2)), 1]
       }
       i<-i+1
    }

  # Neuer Baum
    Tree$Zusammenhangstabelle<- ZHT
    Tree$Varis[,1]           <- 1:length(ZHT_1)
    Tree$CO[,1]              <- 1:length(ZHT_1)

  # Ausgabe:
    return(Tree)
}
