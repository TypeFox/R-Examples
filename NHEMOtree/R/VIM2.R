VIM2 <-
function(Tree, VIM_alt){
  # Vorbereitung
    VIM      <- VIM_alt
    Variablen<- sort(Tree$Varis[,2])

  # Berechnung des VIM fuer einen Baum
    VIM[Variablen-1]<- VIM_alt[Variablen-1]+1/length(Variablen)
    return(VIM)
}
