Mutation_Expansion_NEU <-
function(Tree, Knoten_max_neu, WSK, N_Varis, uS, oS){
  # Zufaellige Auswahl eines Blattes:
    ZHT     <- Tree$Zusammenhangstabelle
    ZHT99   <- which(ZHT[,3]==-99)
    Blatt   <- sample(x=ZHT99, size=1)

  # Maximaler Vaterknoten in ZHT
    MAX<- max(ZHT[,1])

  if (Knoten_max_neu!=0){
  # Erstelle neuen zufaelligen Subtree:
    Wsk_Abnahme<- (WSK-0.001)/Knoten_max_neu
    Subtree    <- Baum_zufaellig_komplett_Sim(Knoten_max=Knoten_max_neu, Wsk=WSK, 
                                              N_Varis=N_Varis, uS=uS, oS=oS)
    ZHT_neu2   <- ZHT_neu      <- Subtree$Zusammenhangstabelle
    Varis_neu  <- Subtree$Varis
    CO_neu     <- Subtree$CO

  # Neuer Subtree benoetigt Eintraege mit anderen Zahlen als Ausgangsbaum
    temp            <- which(ZHT_neu[,3]!=-99)
    ZHT_neu2[,1]    <- ZHT_neu[,1]    +MAX
    ZHT_neu2[temp,3]<- ZHT_neu[temp,3]+MAX
    Varis_neu[,1]   <- Varis_neu[,1]  +MAX
    CO_neu[,1]      <- CO_neu[,1]     +MAX

  # Schreibe neuen Subtree an Stelle des ausgewaehlten Blattes des Ausgangsbaums
    ZHT_mutiert         <- rbind(ZHT, ZHT_neu2)
    ZHT_mutiert[Blatt,3]<- ZHT_neu2[1,1] 
    Varis_mutiert       <- rbind(Tree$Varis, Varis_neu)
    CO_mutiert          <- rbind(Tree$CO,    CO_neu)

  # Ausgabe
    Ausgabe<- list()
    Ausgabe$Zusammenhangstabelle<- ZHT_mutiert
    Ausgabe$Knoten              <- nrow(ZHT_mutiert)/2
    Ausgabe$Varis               <- Varis_mutiert
    Ausgabe$CO                  <- CO_mutiert
    return(Ausgabe)
  }

  # Falls es zu keiner Expansion kommen soll/kann, Ausgabe des Ausgangsbaums!
    if (Knoten_max_neu==0) return(Tree)
}
