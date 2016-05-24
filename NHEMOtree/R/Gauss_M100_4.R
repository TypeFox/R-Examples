Gauss_M100_4 <-
function(Tree, M, gen, alpha, sigma, FKR_max, 
                        Erfolg, Daten, oS, uS, CV.Lauf){
  # Potentieller neuer Baum
    Tree_neu<- Tree

  # Auswahl des Cutoff-Punktes
    Cutoff    <- Tree_neu$CO
    Cutoff_Mut<- 1           # Bei Wurzelbaum, wird der einzige Cutoff mutiert
    if (Tree_neu$Knoten > 1) Cutoff_Mut<- sample(nrow(Cutoff), 1)
    
  # Mutation des Cutoff-Wertes
    Coff<- uS-1

    while(Coff < uS || Coff > oS){ #Solange der Cutoff ausserhalb der Schranken liegt
          Coff<- Cutoff[Cutoff_Mut, 2] + rnorm(n=1, mean=0, sd=sigma)
    }
    Tree_neu$CO[Cutoff_Mut, 2]<- Coff

  # Erfolgreiche Mutation, falls FKR geringer als maximale FKR in Population
    FKR_neu<- 100*FKR_CV(Tree=Tree_neu, Daten=Daten, CV_Laeufe=CV.Lauf)

    if (FKR_neu < FKR_max){
        Erfolg <- Erfolg+1
        Tree   <- Tree_neu   # Bei erfolgreicher Mutation uebernehme diese.
    }

    if (M%%gen == 0) { # Pruefe 1/5-Regel jeweils nach "gen" Generationen
                       # k ist die aktuelle Generation.
      
        # Anpassung von sigma (1/5 Regel)
          if ((Erfolg/gen) > 0.2){
              sigma <- sigma/alpha
          }else{ sigma <- sigma*alpha }

        # Korrektur von sigma (max. sigma=100 erlaubt)
          if (sigma>100) sigma<- 100

        # Zuruecksetzen des Zaehlers    
          Erfolg<- 0  
    }   

  # Anpassen der Baum-Parameter fuer Ausgabe
    Tree$sigma <- sigma
    Tree$Erfolg<- Erfolg

  # Ausgabe:
    return(Tree)
}
