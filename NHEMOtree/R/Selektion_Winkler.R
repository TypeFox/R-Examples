Selektion_Winkler <-
function(Trees){
  N_Parents<- length(Trees)
    
  # Elternteil 1: Roulette-Wheel-Selektion 
    E1<- Roulette_Sim(Trees=Trees)[1]

  # Elternteil 2: zufaellige Selektion
    E2<- E1
    while(E1==E2) E2<- sample(x=1:N_Parents, size=1, replace = FALSE)

  # Ausgabe:
    return(c(E1, E2))
}
