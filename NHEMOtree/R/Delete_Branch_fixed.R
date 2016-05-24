Delete_Branch_fixed <-
function(Tree, Knoten){
#
#    if (length(which(Tree[[1]][,1]==Knoten)) == 0){
#        stop(paste("'Delete_Branch_fixed': Root node does not exist!"))
#    }

  # Teilen der Zusammenhangstabelle in linke und rechte Schritte
    ZHT   <- Tree$Zusammenhangstabelle
    Links <- ZHT[which(ZHT[,2]==0),]
    Rechts<- ZHT[which(ZHT[,2]==1),]

  # Wurzelknoten des zu loeschenden Subtrees
    KC<- c(Knoten, 0)

  # Suche Sohn- und Enkelknoten von KC
    i<- 1; j<-1
    while (j<=length(KC)){
           if (KC[j]!=-99){
           KC[i+1]<- Links[which(Links[,1]==KC[j]),3]
           KC[i+2]<- Rechts[which(Rechts[,1]==KC[j]),3]
           i      <- i+2
           }
           j      <- j+1
    }

  # Alle zu loeschenden Vaterknoten ohne "-99" Eintraege
    KC<- sort(KC)
    while (KC[1]<0) KC<- KC[-1]

  # Tatsaechliches Loeschen aus ZHT, Variablen- und Cutoff-Liste
    temp1<- c(); temp2<- c()
    for (k in 1:length(KC)){
         temp1<- c(temp1, which(ZHT[,1]==KC[k]))
         temp2<- c(temp2, which(Tree$Varis[,1]==KC[k]))
    }
    
    ZHT_kurz  <- ZHT[-temp1,]
    Varis_kurz<- matrix(Tree$Varis[-temp2,], ncol=2)
    CO_kurz   <- matrix(Tree$CO[-temp2,], ncol=2)

  # Loeschen des Tochterknotens mit Wert KC[1] aus ZHT_kurz
    ZHT_kurz[which(ZHT_kurz[,3]==Knoten), 3]<- -99

  # Ausgabe:
    Ausgabe<- list()
    Ausgabe$Zusammenhangstabelle<- ZHT_kurz               # Verkuerzte ZHT
    Ausgabe$Knoten              <- Tree$Knoten-length(KC) # Verringerte Knotenanzahl
    Ausgabe$Varis               <- Varis_kurz             # Verkuerzte Variablentabelle
    Ausgabe$CO                  <- CO_kurz                # Verkuerzte Cutofftabelle
    return(Ausgabe)
}
