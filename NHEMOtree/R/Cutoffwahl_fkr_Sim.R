Cutoffwahl_fkr_Sim <-
function(Tree, Daten, CV.Lauf, AnzCO){
  Tree_neu <- Werte_Reduktion_ZHT2(Tree=Tree)
  ZHT      <- Tree_neu$Zusammenhangstabelle
  Knoten   <- Tree_neu$Knoten
  Varis    <- Tree_neu$Varis
  CO       <- Tree_neu$CO
  CV_Laeufe<- CV.Lauf

  # Auspraegungen
    Auspr      <- sort(unique(ZHT[,1]))
    PCuts      <- calc_N_Pot_Cuts_Sim(Daten=Daten, AnzCO=AnzCO)  # Potentielle Cutoffs aller Variablen
    q          <- list()
    q[[1]]     <- list()
    q[[1]][[1]]<- Auspr[1]
    q[[1]][[2]]<- 1:nrow(Daten)
        
    while (length(q)>0){
      # Erster Knoten ist auf jeden Fall ein Knoten, fuer den ein OptCut bestimmt wird.
      # => Der Status der Kindsknoten kann spaeter abgefragt werden.
        
      # Welcher Knoten aus dem Baum wird jetzt betrachtet?
        vari_index<- q[[1]][[1]]
      
      # Indizes der an diesem Knoten aufzuteilenden Daten (Beobachtungen) 
        akt_index<- q[[1]][[2]]
      
      # Welche Variable aus den Daten gehoert zu dem Knoten im Baum?
        variab<- Varis[which(Varis[,1]==vari_index), 2]

    # Loeschen des aktuellen Knotens aus der Queue
      q[[1]]<- NULL
        
      # Welche potentielle CutOffs hat die aktuelle Variable aus den Daten 
        pc    <- PCuts[[variab]]
      
      # minimale FKR auf "unendlich" setzen
        min_fkr<- 10000000000
        
      # Sampling fuer die Kreuzvalidierung      
      if (CV_Laeufe==nrow(Daten)){ 
        fold<- 1:nrow(Daten)
      }else {
        fold<- sample(x=(1:CV_Laeufe), size=length(akt_index), replace=TRUE)  
      }
      
      # Fuer jeden Poteniellen CutOff wird die FKR bestimmt
      for(i in 1:length(pc)){
        # Initialisierung
          FKR_cv<- c()
          cutoff<- pc[i]
        
        # Fuer jeden Lauf der Kreuzvalidierung
        for(cv in 1:CV_Laeufe){
          
          # Bestimmen der Indexmenge fuer diesen fold
            learning_set_cv<- akt_index[which(fold!=cv)]
            len            <- length(learning_set_cv)
          
            if(len!=0){
               tmp       <- calc_fkr2_Sim(var_index=vari_index, ZHT=ZHT, Varis=Varis, CO=CO,
                                          index_set=learning_set_cv, Daten=Daten, cut=cutoff)
               FKR_cv[cv]<- tmp/len  
            }else{
             FKR_cv[cv]<-0
            }
        }    

      # Mittlere Fehlklassifikationsrate
        fkr<- mean(FKR_cv)

      # Bestimmung des optimalen Cuts
        if(fkr < min_fkr){
          min_fkr<- fkr
          optc   <- cutoff
        }                
      }
      
      # Bestimmung der Indexmenge fuer die naechste Berechnung
        left_index_opt <-intersect(akt_index, which(Daten[, variab]<=optc))
        right_index_opt<-intersect(akt_index, which(Daten[, variab]>optc))
      
      # Setze optimalen Cut
        Tree_neu$CO[which(CO[,1]==vari_index), 2]<- optc
      
      # Hinzufuegen des linken und rechten Knotens
        wurzel<- ZHT[which(ZHT[,1]==vari_index),]
        left  <- wurzel[which(wurzel[,2]==0),]
        right <- wurzel[which(wurzel[,2]==1),]
      
      # Fuege linken Knoten ggf. der Queue hinzu
      if(left[3]!=-99 & length(left_index_opt)!=0){
        q[[(length(q)+1)]]<- list(left[3],left_index_opt)
      } 

      # Fuege rechten Knoten ggf. der Queue hinzu        
      if(right[3]!=-99 & length(right_index_opt)!=0){
        q[[length(q)+1]]<- list(right[3],right_index_opt)
      } 

    } # Ende Schleife fuer alle Knoten

  # Ausgabe:
    return(Tree_neu)
}
