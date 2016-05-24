Fitness_Nachkommen_Sim <-
function(Trees, Daten, Kostenmatrix, CV_Laeufe){
  N_Trees<- length(Trees)

  # Nachkommen 1
    P_N<- Trees[[1]]$Varis[,2]-1

    temp_Na<- sort(P_N); temp_Nb<- temp_Na[1]
    if (length(temp_Na)>1){
        for (k in 1:(length(temp_Na)-1)){
             if (temp_Na[k]!=temp_Na[k+1]) temp_Nb<- c(temp_Nb,temp_Na[k+1])
        }  
    }  
    Rate_N    <- FKR_CV(Tree=Trees[[1]], Daten=Daten, CV_Laeufe=CV_Laeufe)
    Anzahl_N  <- length(temp_Nb)
    Finanzen_N<- getCosts_Tree_sum_Sim(Kostenmatrix=Kostenmatrix, Protein_IDs=temp_Nb)
    FITNESS   <- rbind(Rate_N, Finanzen_N, Anzahl_N)

  # Nachkommen 2 bis N_Trees
    if (N_Trees>1){
    for (i in 2:N_Trees){
         P_N<- Trees[[i]]$Varis[,2]-1
           
         temp_Na<- sort(P_N); temp_Nb<- temp_Na[1]
         if (length(temp_Na)>1){
             for (k in 1:(length(temp_Na)-1)){
                  if (temp_Na[k]!=temp_Na[k+1]) temp_Nb<- c(temp_Nb,temp_Na[k+1])
             }  
         }  
         Rate_N    <- FKR_CV(Tree=Trees[[i]], Daten=Daten, CV_Laeufe=CV_Laeufe)
         Anzahl_N  <- length(temp_Nb)
         Finanzen_N<- getCosts_Tree_sum_Sim(Kostenmatrix=Kostenmatrix, Protein_IDs=temp_Nb)
         Fitness_N <- rbind(Rate_N, Finanzen_N, Anzahl_N)
         FITNESS   <- cbind(FITNESS, Fitness_N)
   }}

   # Ausgabe:
     return(FITNESS)
}
