VIM4_complete <-
function(Tree, Weights, VIM_alt){
  # Vorbereitung
    ZHT      <- Tree$Zusammenhangstabelle
    Varis    <- Tree$Varis
    Varis[,2]<- Varis[,2]-1   # Reduktion um 1, da Spaltenposition der Variablen angegeben ist
    VIM      <- VIM_alt
    Gewichte <- Weights
 
  # Berechnung des VIM fuer einen Baum
    VIM_neu<- VIM4(ZHT, Varis, Gewichte, Knoten=1, Ebene=1, VIM)
    return(VIM_neu)
}
