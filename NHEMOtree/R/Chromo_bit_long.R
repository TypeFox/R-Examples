Chromo_bit_long <-
function(Varis, Gesamtvariablen){
  # Chromosom als Vektor:
    Chromo<- rep(0, length(Gesamtvariablen))
    for (i in 1:length(Varis)) Chromo[which(Gesamtvariablen==Varis[i])] <- 1

  # Chromosom als Zeichenkette
    l_Chromo <- length(Chromo)	
    Chromo2  <- Chromo
    for (j in 1:(l_Chromo-1))  Chromo2[1]<- paste(Chromo2[1], Chromo2[j+1], sep="")
    Chromo3  <- paste(9, Chromo2[1], sep="")

  # Ausgabe als Liste: [[1]] Vektor, [[2]] Zeichenkette mit fuehrender 9
    return(list(Chromo, Chromo3))
}
