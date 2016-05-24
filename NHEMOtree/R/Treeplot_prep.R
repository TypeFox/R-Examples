Treeplot_prep <-
function(Tree, data){
  Tree     <- Werte_Reduktion_ZHT2(Tree)
  Liste    <- list()
  ZHT      <- Tree$Zusammenhangstabelle
  Variablen<- Tree$Varis
  Cutoffs  <- Tree$CO

  # Maximaler Knotenwert in ZHT
    Blatt<- max(ZHT[,1])

  # Pro Vaterknoten eine Liste
    for(i in 1:(nrow(ZHT)/2)){
        j<- 2*i-1 # alle ungeraden Zeilen

        # Beide Kinder keine Blaetter
        if (ZHT[j,3]!=-99 && ZHT[j+1,3]!=-99){ 
            Liste[[i]]<- list(id = as.integer(ZHT[j,1]), 
                           split = partysplit(varid = as.integer(Variablen[which(Variablen[,1]==ZHT[j,1]),2]), breaks = Cutoffs[which(Cutoffs[,1]==ZHT[j,1]), 2]), 
                            kids = c(as.integer(ZHT[j,3]), as.integer(ZHT[j+1,3])))
        }

        # Linkes Kind ist Blatt, rechtes Kind kein Blatt
        if (ZHT[j,3]==-99 && ZHT[j+1,3]!=-99){ # Linkes Kind ist Blatt, rechtes Kind kein Blatt
            Liste[[i]]<- list(id = as.integer(ZHT[j,1]), 
                           split = partysplit(varid = as.integer(Variablen[which(Variablen[,1]==ZHT[j,1]), 2]), breaks = Cutoffs[which(Cutoffs[,1]==ZHT[j,1]), 2]), 
                            kids = c(as.integer(max(Blatt)+1), as.integer(ZHT[j+1,3])))
            Blatt     <- c(Blatt, max(Blatt)+1)
        }

        # Linkes Kind kein Blatt, rechtes Kind ist Blatt
        if (ZHT[j,3]!=-99 && ZHT[j+1,3]==-99){ 
            Liste[[i]]<- list(id = as.integer(ZHT[j,1]), 
                           split = partysplit(varid = as.integer(Variablen[which(Variablen[,1]==ZHT[j,1]), 2]), breaks = Cutoffs[which(Cutoffs[,1]==ZHT[j,1]), 2]), 
                            kids = c(as.integer(ZHT[j,3]), as.integer(max(Blatt)+1)))
            Blatt     <- c(Blatt, max(Blatt)+1)
        }

        # Beide Kinder Blaetter
        if (ZHT[j,3]==-99 && ZHT[j+1,3]==-99){ # Beide Kinder Blaetter
            Liste[[i]]<- list(id = as.integer(ZHT[j,1]), 
                           split = partysplit(varid = as.integer(Variablen[which(Variablen[,1]==ZHT[j,1]), 2]), breaks = Cutoffs[which(Cutoffs[,1]==ZHT[j,1]), 2]), 
                            kids = c(as.integer(max(Blatt)+1), as.integer(max(Blatt)+2)))
            Blatt     <- c(Blatt, max(Blatt)+1, max(Blatt)+2)
        }
    }#Liste aller Vaterknoten fertig!

  # Pro Blatt eine Liste:
    Blatt<- Blatt[-1]
    for (l in Blatt) Liste[[l]]<- list(id = l)

  # Erstellen des Baumes
    Knotenliste<- as.partynode(Liste)
    Baum_malen <- party(Knotenliste, data=data) 

  # Ausgabe:
    return(Baum_malen)
}
