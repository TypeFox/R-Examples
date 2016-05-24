Vater_in_ZHT_prozentual <-
function(ZHT, zaehler, Knoten_max, Wsk, Wsk_Abnahme){
  Vaterknoten<- 1
  while(zaehler>Vaterknoten && Vaterknoten<=Knoten_max){
    Vaterknoten<- Vaterknoten+1

    # Links vom Vaterknoten
      if (runif(1, min=0, max=1)<= Wsk){
          zaehler<- zaehler+1
          Z      <- c(Vaterknoten, 0, zaehler)
          ZHT    <- rbind(ZHT, Z)
      } else{Z   <- c(Vaterknoten, 0, -99)
             ZHT <- rbind(ZHT, Z)}

    # Rechts vom Vaterknoten
      if (runif(1, min=0, max=1)<= Wsk){
          zaehler<- zaehler+1
          Z      <- c(Vaterknoten, 1, zaehler)
          ZHT    <- rbind(ZHT, Z)
      } else{Z   <- c(Vaterknoten, 1, -99)
             ZHT <- rbind(ZHT, Z)}

    Wsk<- Wsk-Wsk*Wsk_Abnahme
  }

  # Falls letzter Tochterknoten>Knoten_max, ersetze durch -99
    while(max(ZHT[,3])>Knoten_max) ZHT[which.max(ZHT[,3]), 3]<- -99

  # Falls letzter Vaterknoten unvollstaendig, d.h. kein rechts, dann einfuegen:
    if (ZHT[nrow(ZHT),2]==0){
        Z_1<- c(max(ZHT[,1]), 1, -99)
        ZHT<- rbind(ZHT, Z_1)
    }

  # Falls letzter Vaterknoten kleiner als letzter Tochterknoten, fuelle ZHT auf
    if (max(ZHT[,1]) < max(ZHT[,3])){
        temp<- (max(ZHT[,1])+1):max(ZHT[,3])
        Z1  <- sort(c(temp,temp))
        Z2  <- rep(c(0,1), length(temp))
        Z3  <- rep(-99, 2*length(temp))
        Z   <- cbind(Z1, Z2, Z3)
        
        ZHT<- rbind(ZHT, Z)
    }

  # Loeschen, falls zu viele Vaterknoten in ZHT sind
    while(ZHT[nrow(ZHT),] > Knoten_max && ZHT[nrow(ZHT), 3]==-99) ZHT<- ZHT[-nrow(ZHT),]
    
  # Ausgabe:        
    Ausgabe                     <- list()
    Ausgabe$Zusammenhangstabelle<- ZHT 
    Ausgabe$Knoten              <- max(ZHT[,1])
    return(Ausgabe)
}
