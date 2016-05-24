FKR_CV <-
function(Tree, Daten, CV_Laeufe){

  Tree     <- Werte_Reduktion_ZHT2(Tree=Tree)  # Korrekte Nummerierung der ZHT-Spalteneintraege
  ZHT      <- Tree$Zusammenhangstabelle
  N_Knoten <- Tree$Knoten
  Variablen<- Tree$Varis
  Cutoffs  <- Tree$CO
  BeobNeu  <- list()
  FKR.cv   <- c()

  ##########################
  # Keine Kreuzvalidierung #
  ##########################
    if (CV_Laeufe==1){
        FKR<- FKR_fast_NEU2(Tree=Tree, Daten=Daten)
        return(FKR)
    }

  ####################
  # Kreuzvalidierung #
  ####################
    fold<- sample(x=(1:CV_Laeufe), size=nrow(Daten), replace=TRUE)
    if (CV_Laeufe==nrow(Daten)) fold<- 1:nrow(Daten)

    for(cv in 1:CV_Laeufe){
        Learning_set.cv<- Daten[which(fold!=cv),]

        FKR.cv[cv]<- FKR_fast_NEU2(Tree=Tree, Daten=Learning_set.cv)  # Fehlklassifikationsrate
    }

    # Ausgabe:
      return(mean(FKR.cv))
}
