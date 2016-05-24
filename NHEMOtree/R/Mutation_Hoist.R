Mutation_Hoist <-
function(Tree, KC=0){
  # Falls Tree Wurzelbaum, keine Hoist-Mutation moeglich:
    if (Tree$Knoten==1) return(Tree)

  # Teilen der Zusammenhangstabelle in linke und rechte Schritte
    ZHT   <- Tree$Zusammenhangstabelle
    Links <- ZHT[which(ZHT[,2]==0),]
    Rechts<- ZHT[which(ZHT[,2]==1),]

  # Zufaellige Auswahl eines Vaterknotes des Subbaums, 
  # der neuer Baum werden soll (ausser Wurzelknoten)
    if (KC==0){
        KC_temp<- sample(2:Tree$Knoten, 1)
        KC     <- ZHT[2*KC_temp, 1]
    }
    KC[2]<- 0

  # Suche Sohn- und Enkelknoten des Subbaums
    i<- 1; j<-1
    while (j<=length(KC)){
      if (KC[j]!=-99){
          KC[i+1]<- Links[which(Links[,1]==KC[j]),3]
          KC[i+2]<- Rechts[which(Rechts[,1]==KC[j]),3]
          i      <- i+2
      }
      j<- j+1
    }

  # Alle zu loeschenden Vaterknoten ohne "-99" Eintraege
    KC<- sort(KC)
    while (KC[1]<0) KC<- KC[-1]

  # Loeschen aller Knoten etc aus ZHT, Variablen- und Cutofftabelle,
  # die nicht in Subbaum enthalten sind
    temp1<- c(); temp2<- c()
    for (k in 1:length(KC)){
         temp1<- c(temp1, which(ZHT[,1]==KC[k]))
         temp2<- c(temp2, which(Tree$Varis[,1]==KC[k]))
    }

    ZHT_Subbaum  <- ZHT[temp1,]
    Varis_Subbaum<- matrix(Tree$Varis[temp2,], ncol=2)
    CO_Subbaum   <- matrix(Tree$CO[temp2,], ncol=2)

  # Ausgabe:
    Ausgabe                     <- list()
    Ausgabe$Zusammenhangstabelle<- ZHT_Subbaum
    Ausgabe$Knoten              <- nrow(ZHT_Subbaum)/2
    Ausgabe$Varis               <- Varis_Subbaum    
    Ausgabe$CO                  <- CO_Subbaum
    Ausgabe$geloescht           <- KC[1]
    Ausgabe                     <- Werte_Reduktion_ZHT2(Tree=Ausgabe)
    return(Ausgabe)
}
