"checkWW" <- 
function(X,WW = NULL,ngroups = NULL,ord=FALSE) {

# Diese Funktion ueberprueft einen Vektor der Liste "indices" auf vorhandene Wechselwirkungen (WW)
# und sortiert dabei Permutationen aus, indem WW vor dem Hauptwirkungen positioniert sind.
# Benoetigt wird eine Matrix mit 2 Spalten, wobei in der ersten Spalte die Position der WW im originaldatensatz
# angegeben ist und in der 2. Spalte die Position einer benoetigten Hauptkomponenten angegeben ist.
# Die Matrix WW enthaelt somit alle benoetigten Hauptkompunenten zu jeder Wechselwirkung   

## X ist Vektor

    erg <- TRUE
    if (!is.null(WW)){
      if (!is.null(ngroups)) 
           X <- list2vec(ngroups[X])-1  ## effect numbers (intercept excluded)
        for (i in 1:nrow(WW)) {
            # passt im Fall von Gruppen die Reihenfolge an
            # Prueft, ob eine WW vorhanden ist.
            if (any(X == WW[i,1]) && erg){
                if (!WW[i,2] %in% X) erg <- FALSE
                # Prueft, ob die zugehoerige HK vorhanden ist.
                if (ord){
                    # Prueft, ob die WW nach der HK positioniert ist.
                    if (which(X == WW[i,1]) < which(X == WW[i,2]))
                        erg <- FALSE
                        }
            }
        }

      }

    return(erg)
}
