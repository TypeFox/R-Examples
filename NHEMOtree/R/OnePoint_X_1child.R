OnePoint_X_1child <-
function(Eltern_Chromos, crossover_prob){
# if (nrow(Eltern_Chromos)!=2) stop(paste("Falsche Dimension des Inputs in 'OnePoint_X_1child'"))
  E1<- Eltern_Chromos[1,]
  E2<- Eltern_Chromos[2,]
  N <- length(E1)

  temp<- runif(n=1, min=-0.01, max=1.01)

  # Falls Rekombination stattfindet:
    if (temp <= crossover_prob){
      Position<- round(runif(n=1, min=0.5, max=(N-0.51)))

    # Erstellung der Nachkommen aus den Eltern
      C1_1<- E1[1:Position]
      C1_2<- E2[(Position+1):N]
      C1  <- c(C1_1, C1_2)

      C2_1<- E2[1:Position]
      C2_2<- E1[(Position+1):N]
      C2  <- c(C2_1, C2_2)

    # Welcher Nachkommen soll ueberleben?
      temp2<- runif(n=1, min=0, max=1)
              if (temp2 <= 0.5) Child<- C1
              if (temp2 >  0.5) Child<- C2
    } # Falls Zufallszahl kleiner als Rekombinationsrate --> Rekombination!

  # Falls KEINE Rekombination stattfindet:
    if(temp > crossover_prob){
      # Welcher Nachkommen soll ueberleben?
        temp2<- runif(n=1, min=0, max=1)
                if (temp2 <= 0.5) Child<- E1
                if (temp2 >  0.5) Child<- E2

    } # Falls Zufallszahl groesser als Rekombinationsrate --> Keine Rekombination!!!

  # Ausgabe:
    return(Child)
}
