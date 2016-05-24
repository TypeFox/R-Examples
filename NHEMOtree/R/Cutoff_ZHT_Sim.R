Cutoff_ZHT_Sim <-
function(Knoten, uS, oS){
  CO <- rnorm(n=Knoten, mean=(oS-uS)/2, sd=(oS-uS)/2*0.1)
  # Ausgabe:
    return(cbind(1:Knoten, CO))
}
