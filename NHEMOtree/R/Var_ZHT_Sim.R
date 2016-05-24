Var_ZHT_Sim <-
function(Knoten, N_Varis){
  K      <- c() 
  for (i in 1:Knoten) K[i] <- sample(2:(N_Varis+1), 1)  # Variablenposition in R-Tabelle
  return(cbind(1:Knoten, K))
}
