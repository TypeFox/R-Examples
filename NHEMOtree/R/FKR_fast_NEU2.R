FKR_fast_NEU2 <-
function(Tree, Daten){

  Tree     <- Werte_Reduktion_ZHT2(Tree)  # Korrekte Nummerierung der ZHT-Spalteneintraege
  ZHT      <- Tree$Zusammenhangstabelle
  N_Knoten <- Tree$Knoten
  Variablen<- Tree$Varis
  Cutoffs  <- Tree$CO
  BeobNeu  <- list()
  FK       <- c()                         # Fehlklassifikation pro Analysepunkt (absolute Zahlen)

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
            r<- ZHT[which(ZHT[,3]==i),]        # Zeilen mit 3.Spalte = i 
            k<- 2*r[1]-(1-r[2])               
            BeobNeu[[2*i-1]]<- intersect(which(Daten[, Variablen[i,2]] <= Cutoffs[i,2]), BeobNeu[[k]]) 
            BeobNeu[[2*i]]  <- intersect(which(Daten[, Variablen[i,2]] >  Cutoffs[i,2]), BeobNeu[[k]])
        }
     }

  # Auswertung der Fehlklassifikationsrate
    T<- ZHT[which(ZHT[,3]==-99),]   # Alle Blaetter

    l<- 1
    while(l <= nrow(T)){
          r<- T[l,]
          k<- 2*r[1]-(1-r[2])
          SCHNITT<- BeobNeu[[k]]
          Subtyp <- Daten[SCHNITT,1]   
          FK[l]  <- length(Subtyp)-max(table(Subtyp))
          l      <- l+1
    }

    FKR<- sum(FK)/nrow(Daten)  # misclassification rate

  # Ausgabe:
    return(FKR)
}
