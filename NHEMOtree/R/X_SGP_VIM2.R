X_SGP_VIM2 <-
function(E1, E2, N_Varis, uS, oS, 
                      Max_Knoten, Daten, CV_Laeufe, VIM){

  # Falls beide Baeume Wurzelbaeume, direkte Rueckgabe
    if (E1$Knoten==1 || E2$Knoten==1) return(list(E1,E2))

  # Zunaechst Betrachtung des groesseren Baums
    if (E2$Knoten > E1$Knoten){
        temp<- E1
        E1  <- E2
        E2  <- temp
    }

  # Abbruch, falls Knotenanzahl beider Baeume > 2 x Max_Knoten
    if (E1$Knoten+E2$Knoten > 2*Max_Knoten) stop(paste("'X_SGP_VIM2': The trees are too big!"))

  ZHT1<- E1$Zusammenhangstabelle
  ZHT2<- E2$Zusammenhangstabelle

  #################
  # 1. Elternbaum #
  #################
  # Auswahl des Rekombinationspunktes des 1. Elternteils (ausser Wurzelknoten):
    KC1   <- WK_Auswahl_Rang(Tree=E1, VIM=VIM, N_Varis=N_Varis)
    KC1[2]<- 0

  # Suche Sohn- und Enkelknoten der Subtrees
    Links_1 <- ZHT1[which(ZHT1[,2]==0),]
    Rechts_1<- ZHT1[which(ZHT1[,2]==1),]
    i<- 1; j<-1
    while (j<=length(KC1)){
      if (KC1[j]!=-99){
          KC1[i+1]<- Links_1[which(Links_1[,1]==KC1[j]),3]
          KC1[i+2]<- Rechts_1[which(Rechts_1[,1]==KC1[j]),3]
          i       <- i+2
      }
      j<- j+1
    }
 
  # Alle zu loeschenden Vaterknoten ohne "-99" Eintraege
    KC1<- sort(KC1)
    while (KC1[1]<0) KC1<- KC1[-1]

  # Extrahieren der Subtrees und Loeschen der Subtrees aus den Eltern-ZHT:
    temp1a<- c(); temp1b<- c()
    for (i in 1:length(KC1)){
         temp1a<- c(temp1a, which(ZHT1[,1]==KC1[i]))     # Zu loeschende Zeilen
         temp1b<- c(temp1b, which(E1$Varis[,1]==KC1[i])) # Zu loeschende Zeilen
    }

    ZHT1_kurz          <- ZHT1[-temp1a,]
    Varis1_kurz        <- matrix(E1$Varis[-temp1b,], ncol=2)
    CO1_kurz           <- matrix(E1$CO[-temp1b,], ncol=2)

    ZHT1_Subtree_sauber<- ZHT1_Subtree<- ZHT1[temp1a,]
    Varis1_Subtree     <- matrix(E1$Varis[temp1b,], ncol=2)
    CO1_Subtree        <- matrix(E1$CO[temp1b,], ncol=2)


  #################
  # 2. Elternbaum #
  #################
  # Auswahl des Rekombinationspunktes des 2. Elternteils (ausser Wurzelknoten):
    KC2     <- WK_Auswahl_Rang(Tree=E2, VIM=VIM, N_Varis=N_Varis)
    KC2[2]<- 0

  # Suche Sohn- und Enkelknoten der Subtrees
    Links_2 <- ZHT2[which(ZHT2[,2]==0),]
    Rechts_2<- ZHT2[which(ZHT2[,2]==1),]
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

  # Extrahieren der Subtrees und Loeschen der Subtrees aus den Eltern-ZHT
    temp2a<- c(); temp2b<- c()
    for (i in 1:length(KC2)){
         temp2a<- c(temp2a, which(ZHT2[,1]==KC2[i]))
         temp2b<- c(temp2b, which(E2$Varis[,1]==KC2[i]))
    }

    ZHT2_kurz          <- ZHT2[-temp2a,]
    Varis2_kurz        <- matrix(E2$Varis[-temp2b,], ncol=2)
    CO2_kurz           <- matrix(E2$CO[-temp2b,], ncol=2)
 
    ZHT2_Subtree_sauber<- ZHT2_Subtree<- ZHT2[temp2a,]
    Varis2_Subtree     <- matrix(E2$Varis[temp2b,], ncol=2)
    CO2_Subtree        <- matrix(E2$CO[temp2b,], ncol=2)


  ##########################
  # KOMBINATION DER ELTERN #
  ##########################
  # Elternteil 1 + Subtree von Elternteil 2
    MAX1<- max(ZHT1_kurz[,1])  # Maximaler Vaterknoten in der verkuerzten ZHT1
    
    # Subtree2 benoetigt Eintraege mit anderen Zahlen als ZHT1_kurz
      temp99_2                       <- which(ZHT2_Subtree_sauber[,3]!=-99)
      ZHT2_Subtree_sauber[,1]        <- ZHT2_Subtree_sauber[,1]        -KC2[1]+1+MAX1
      ZHT2_Subtree_sauber[temp99_2,3]<- ZHT2_Subtree_sauber[temp99_2,3]-KC2[1]+1+MAX1
      Varis2_Subtree[,1]             <- Varis2_Subtree[,1]             -KC2[1]+1+MAX1
      CO2_Subtree[,1]                <- CO2_Subtree[,1]                -KC2[1]+1+MAX1

    # Schreibe neuen Subtree an Stelle des alten Subtrees: 
      ZHT1_kurz[which(ZHT1_kurz[,3]==KC1[1]),3]<- ZHT2_Subtree_sauber[1,1] 
      ZHT_child1                               <- rbind(ZHT1_kurz, ZHT2_Subtree_sauber)
      Varis_child1                             <- rbind(Varis1_kurz, Varis2_Subtree)
      CO_child1                                <- rbind(CO1_kurz, CO2_Subtree)
      
      Child1<- list(Zusammenhangstabelle=ZHT_child1, Knoten=nrow(ZHT_child1)/2, 
                    Varis=Varis_child1, CO=CO_child1)

  # Elternteil 2 + Subtree von Elternteil 1
    MAX2<- max(ZHT2_kurz[,1])  # Maximaler Vaterknoten in der verkuerzten ZHT2
    
    # Subtree2 benoetigt Eintraege mit anderen Zahlen als ZHT2_kurz
      temp99_1                       <- which(ZHT1_Subtree_sauber[,3]!=-99)
      ZHT1_Subtree_sauber[,1]        <- ZHT1_Subtree_sauber[,1]        -KC1[1]+1+MAX2
      ZHT1_Subtree_sauber[temp99_1,3]<- ZHT1_Subtree_sauber[temp99_1,3]-KC1[1]+1+MAX2
      Varis1_Subtree[,1]             <- Varis1_Subtree[,1]             -KC1[1]+1+MAX2
      CO1_Subtree[,1]                <- CO1_Subtree[,1]                -KC1[1]+1+MAX2

    # Schreibe neuen Subtree an Stelle des alten Subtrees: 
      ZHT2_kurz[which(ZHT2_kurz[,3]==KC2[1]),3]<- ZHT1_Subtree_sauber[1,1] 
      ZHT_child2                               <- rbind(ZHT2_kurz, ZHT1_Subtree_sauber)
      Varis_child2                             <- rbind(Varis2_kurz, Varis1_Subtree)
      CO_child2                                <- rbind(CO2_kurz, CO1_Subtree)
      
      Child2<- list(Zusammenhangstabelle=ZHT_child2, Knoten=nrow(ZHT_child2)/2, 
                    Varis=Varis_child2, CO=CO_child2)

  # Fitnessorientierte Hoist-Mutation, falls Kinder groesser als maximal erlaubt
    if (Child1$Knoten>Max_Knoten) Child1<- Delete_EL(Tree=Child1, Daten=Daten, N_Varis=N_Varis, uS=uS, oS=oS)
    if (Child1$Knoten>Max_Knoten){
        Child1<- Mutation_Hoist_Fitness(Tree=Child1, Datensatz=Daten, CV_Laeufe=CV_Laeufe)
        while (Child1$Knoten>Max_Knoten) Child1<- Delete_Branch_random(Tree=Child1)  
    }

    if (Child2$Knoten>Max_Knoten) Child2<- Delete_EL(Tree=Child2, Daten=Daten, N_Varis=N_Varis, uS=uS, oS=oS)
    if (Child2$Knoten>Max_Knoten){
        Child2<- Mutation_Hoist_Fitness(Tree=Child2, Datensatz=Daten, CV_Laeufe=CV_Laeufe)
        while (Child2$Knoten>Max_Knoten) Child2<- Delete_Branch_random(Tree=Child2)  
    }

  # Ausgabe:
    Ausgabe     <- list()
    Ausgabe[[1]]<- Child1
    Ausgabe[[2]]<- Child2
    return(Ausgabe)
}
