Mutation_Branch_fully_Sim <-
function(Tree, Knoten_max, WSK, Wsk_Abnahme, N_Varis, uS, oS){
 
  # Korrekte Nummerierung der ZHT-Spalteneintraege
    Tree<- Werte_Reduktion_ZHT2(Tree=Tree) 

  ########################
  # Loeschen eines Astes #
  ########################
  # Falls Tree Wurzelbaum, nur Auswahl der Seite zum Hinzufuegen moeglich
  if (Tree$Knoten==1){
      ZHT_kurz                     <- Tree$Zusammenhangstabelle
      ZHT_kurz[sample(c(1,2), 1),3]<- 2
      KC                           <- 2
      Varis_kurz                   <- Tree$Varis
      CO_kurz                      <- Tree$CO
  }

  # Falls Tree kein Wurzelbaum:
  if (Tree$Knoten>1){ 
  # Teilen der Zusammenhangstabelle in linke und rechte Schritte
    ZHT   <- Tree$Zusammenhangstabelle
    Links <- ZHT[which(ZHT[,2]==0),]
    Rechts<- ZHT[which(ZHT[,2]==1),]

  # Zufaellige Auswahl eines Vaterknotes des zu loeschenden Astes 
  # (ausser Wurzelknoten des Originalbaumes)
    KC_temp<- sample(2:Tree$Knoten, 1)
    if (Tree$Knoten==2) KC_temp<- 2
    KC     <- ZHT[2*KC_temp, 1]
    KC[2]  <- 0

  # Suche Sohn- und Enkelknoten von KC
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

  # Tatsaechliches Loeschen aus ZHT, Variablen- und Cutofftabelle
    temp1<- c(); temp2<- c()
    for (k in 1:length(KC)){
         temp1<- c(temp1, which(ZHT[,1]==KC[k]))
         temp2<- c(temp2, which(Tree$Varis[,1]==KC[k]))
    }
    ZHT_kurz  <- ZHT[-temp1,]
    Varis_kurz<- matrix(Tree$Varis[-temp2,], ncol=2)
    CO_kurz   <- matrix(Tree$CO[-temp2,], ncol=2)
  }

  ###############################
  # Erstellen eines neuen Astes #
  ###############################
  Knoten_kurz<- nrow(ZHT_kurz)/2        # Groesse des gestutzten Baumes
  Knoten_neu <- Knoten_max-Knoten_kurz  # Maximales Groesse des Subbaumes

  if (Knoten_neu<=0){
      ZHT_kurz[which(ZHT_kurz[,3]==KC[1]),3]<- -99
 
      Ausgabe<- list()
      Ausgabe$Zusammenhangstabelle<- ZHT_kurz     # Verkuerzte ZHT
      Ausgabe$Knoten              <- Knoten_kurz  # Verringerte Knotenanzahl
      Ausgabe$Varis               <- Varis_kurz   # Verkuerzte Variablentabelle
      Ausgabe$CO                  <- CO_kurz      # Verkuerzte Cutofftabelle
      Ausgabe$Loeschposition      <- KC[1]        # Vaterknotes des geloeschten Astes
      return(Ausgabe)
  }

  if (Knoten_neu>0){
  # Maximaler Wert eines Vaterknotens in verkuerzter ZHT
    MAX<- max(ZHT_kurz[,1])

  # Erstelle neuen zufaelligen Subtree
    Subtree  <- Baum_zufaellig_komplett_Sim(Knoten_max=Knoten_neu, Wsk=WSK, 
                                            N_Varis=N_Varis, uS=uS, oS=oS)
    ZHT_neu2 <- ZHT_neu  <- Subtree$Zusammenhangstabelle
    Varis_neu<- Subtree$Varis
    CO_neu   <- Subtree$CO

    temp            <- which(ZHT_neu[,3]!=-99)
    ZHT_neu2[,1]    <- ZHT_neu[,1]    +MAX
    ZHT_neu2[temp,3]<- ZHT_neu[temp,3]+MAX
    Varis_neu[,1]   <- Varis_neu[,1]  +MAX
    CO_neu[,1]      <- CO_neu[,1]     +MAX

  # Schreibe neuen Subtree an Stelle des alten Subtrees
    ZHT_mutiert  <- rbind(ZHT_kurz, ZHT_neu2)
    ZHT_mutiert[which(ZHT_mutiert[,3]==KC[1]),3]<- ZHT_neu2[1,1] # Wurzelknoten von ZHT_neu, wo alter Vaterknoten als Tochterknoten
    Varis_mutiert<- rbind(Varis_kurz, Varis_neu)
    CO_mutiert   <- rbind(CO_kurz,    CO_neu)

  # Ausgabe:
    Ausgabe                     <- list()
    Ausgabe$Zusammenhangstabelle<- ZHT_mutiert          # Verkuerzte ZHT
    Ausgabe$Knoten              <- nrow(ZHT_mutiert)/2  # Verringerte Knotenanzahl
    Ausgabe$Varis               <- Varis_mutiert        # Verkuerzte Variablentabelle
    Ausgabe$CO                  <- CO_mutiert           # Verkuerzte Cutofftabelle
    Ausgabe$Loeschposition      <- KC[1]                # Vaterknotes des geloeschten Astes
    return(Ausgabe)
  }
}
