X_1Point_GP_WK <-
function(E1, E2){

  # Wurzelbaum als erster Elternbaum
    if (E2[[2]]==1){
        TEMP_E1<- E2
        E2     <- E1
        E1     <- TEMP_E1
    }
    ZHT1<- E1$Zusammenhangstabelle
    ZHT2<- E2$Zusammenhangstabelle

  ##############################
  # 1. Elternbaum - Wurzelbaum #
  ##############################
    KC1     <- sample(1:2, 1) # Erweiterung des Wurzelbaums links oder rechts

  ###################################
  # 2. Elternbaum - Kein Wurzelbaum #
  ###################################
  # Teilen der Zusammenhangstabelle in linke und rechte Schritte
    Links_2 <- ZHT2[which(ZHT2[,2]==0),]
    Rechts_2<- ZHT2[which(ZHT2[,2]==1),]

  # Auswahl des Subbaums links oder rechts des Wurzelknoten
    RP2<- sort(ZHT2[1:2,3])           # (*)
    KC2<- sample(RP2, 1); KC2[2]<- 0  # Rekombinationsknoten des 2. Elternteils

    # Falls Baum nur zu einer Seite waechst und -99 ausgewaehlt wird
    # --> KEINE REKOMBINATION
    if (KC2[1]==-99){
        Child1<- E1; Child1<- Child1[-5]  # Loeschen der elterlichen Lauf-ID
        Child2<- E2; Child2<- Child2[-5]  # Loeschen der elterlichen Lauf-ID
    }

  if (KC2[1]!=-99){
  # Suche Sohn- und Enkelknoten der Subtrees
    i<- 1; j<-1
    while (j<=length(KC2)){
           if (KC2[j]!=-99){
           KC2[i+1]<- Links_2[which(Links_2[,1]==KC2[j]),3]
           KC2[i+2]<- Rechts_2[which(Rechts_2[,1]==KC2[j]),3]
           i       <- i+2
           }
           j<- j+1
    }

  # Alle zu loeschenden Vaterknoten ohne "-99" Eintraege
    KC2<- sort(KC2)
    while (KC2[1]<0) KC2<- KC2[-1]

  # Extrahieren der Subtrees und Loeschen der Subtrees aus den Eltern-ZHT:
    tempa<- c(); tempb<- c()
    for (i in 1:length(KC2)){
         tempa<- c(tempa, which(ZHT2[,1]==KC2[i]))
         tempb<- c(tempb, which(E2$Varis[,1]==KC2[i]))
    }

    ZHT2_kurz          <- ZHT2[-tempa,]
    Varis2_kurz        <- matrix(E2$Varis[-tempb,], ncol=2)
    CO2_kurz           <- matrix(E2$CO[-tempb,], ncol=2)
    ZHT2_Subtree_sauber<- ZHT2_Subtree<- ZHT2[tempa,]
    Varis2_Subtree     <- matrix(E2$Varis[tempb,], ncol=2)
    CO2_Subtree        <- matrix(E2$CO[tempb,], ncol=2)


  ##########################
  # KOMBINATION DER ELTERN #
  ##########################
  # Wurzelbaum + Subtree von Elternteil 2
    # Subtree2 benoetigt Eintraege mit anderen Zahlen als ZHT1_kurz
      temp99_2                       <- which(ZHT2_Subtree_sauber[,3]!=-99)
      ZHT2_Subtree_sauber[,1]        <- ZHT2_Subtree_sauber[,1]-KC2[1]        +2 
      ZHT2_Subtree_sauber[temp99_2,3]<- ZHT2_Subtree_sauber[temp99_2,3]-KC2[1]+2 
      Varis2_Subtree[,1]             <- Varis2_Subtree[,1]             -KC2[1]+2
      CO2_Subtree[,1]                <- CO2_Subtree[,1]                -KC2[1]+2
 
    # Schreibe neuen Subtree an Stelle des Wurzelbaums: 
      ZHT1[KC1,3] <- ZHT2_Subtree_sauber[1,1] 
      ZHT_child1  <- rbind(ZHT1, ZHT2_Subtree_sauber)
      Varis_child1<- rbind(E1$Varis, Varis2_Subtree)
      CO_child1   <- rbind(E1$CO, CO2_Subtree)
      
      Child1<- list(Zusammenhangstabelle=ZHT_child1, Knoten=nrow(ZHT_child1)/2, 
                    Varis=Varis_child1, CO=CO_child1)

  # Elternteil 2 an die Stelle der Subtree-Abtrennung eine '-99'
    ZHT_child2<- ZHT2_kurz
    ZHT_child2[which(ZHT_child2[,3]==KC2[1]),3]<- -99

    Child2<- list(Zusammenhangstabelle=ZHT_child2, Knoten=nrow(ZHT_child2)/2, 
                  Varis=Varis2_kurz, CO=CO2_kurz)
  } # if (KC2[1]!=-99)
   
  # Ausgabe:
    Ausgabe<- list()
    Ausgabe[[1]]<- Child1
    Ausgabe[[2]]<- Child2
    return(Ausgabe)
}
