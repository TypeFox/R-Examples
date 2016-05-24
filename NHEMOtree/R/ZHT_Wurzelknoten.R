ZHT_Wurzelknoten <-
function(Wsk){
  ZHT    <- matrix(c(1,0,0), nrow=1, byrow=T)
  zaehler<- 1

  # Links vom Wurzelknoten
    if (runif(1, min=0, max=1)<= Wsk){
        zaehler  <- zaehler+1
        ZHT[1, 3]<- zaehler
    } else ZHT[1,3]<- -99

  # Rechts vom Vaterknoten
    if (runif(1, min=0, max=1)<= Wsk){
        zaehler<- zaehler+1
        Z      <- c(1, 1, zaehler)
        ZHT    <- rbind(ZHT, Z)
    } else{Z<- c(1, 1, -99)
           ZHT<- rbind(ZHT, Z)}

  Ausgabe<- list()
  Ausgabe$Zusammenhangstabelle<-ZHT  # Ausgabe der Zusammenhangstabelle
  Ausgabe$Knoten<-zaehler            # Anzahl der Knoten im Baum
  return(Ausgabe)
}
