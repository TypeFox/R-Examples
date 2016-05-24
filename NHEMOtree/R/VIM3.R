VIM3 <-
function(Tree, VIM_alt){
  # Vorbereitung
    VIM      <- VIM_alt
    Variablen<- sort(Tree$Varis[,2])
    Anzahl   <- length(Variablen)

  # Berechnung des VIM fuer einen Baum
    while(length(Variablen)>0){
          Variablen2<- c()
          if (length(Variablen)>1){
              for (l in 1:(length(Variablen)-1)){
                   if (Variablen[l]==Variablen[l+1]) Variablen2<- c(Variablen2, Variablen[l+1])
          }}  
          VIM[Variablen-1]<- VIM[Variablen-1]+1/Anzahl
          Variablen       <- Variablen2
    }
    return(VIM)
}
