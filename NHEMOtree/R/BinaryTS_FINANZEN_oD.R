BinaryTS_FINANZEN_oD <-
function(N_Parents, Fitness, Vectors){
  P1<- P2<- 0
  while (P1==P2){
         Eltern1<- round(runif(n=2, min=0.5, max=N_Parents+0.49))
         Eltern2<- round(runif(n=2, min=0.5, max=N_Parents+0.49))
 
       # 1. Elternteil:
         Elternpaar1<- Fitness[1:2, Eltern1]
         Vergleich  <- nds_rank(Elternpaar1)

         if (Vergleich[1]<=Vergleich[2]) P1 <- Eltern1[1]
         if (Vergleich[1]>Vergleich[2])  P1 <- Eltern1[2]

       # 2. Elternteil:
         Elternpaar2<- Fitness[1:2, Eltern2]
         Vergleich  <- nds_rank(Elternpaar2)

         if (Vergleich[1]<=Vergleich[2]) P2 <- Eltern2[1]
         if (Vergleich[1]>Vergleich[2])  P2 <- Eltern2[2]
       }
       P<- c(P1, P2)

  # Ausgabe der Elternchromosomen und deren Position in 'Ind'
    return(list(Vectors[P,], P))
}
