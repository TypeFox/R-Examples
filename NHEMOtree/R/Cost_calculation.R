Cost_calculation <-
function(Varis, KOSTEN){

  Varis    <- unique(Varis)
  Varis_num<- c()
  for (i in 1: length(Varis)) Varis_num[i]<- which(KOSTEN[,1]==Varis[i])
  
  Total_costs<- sum(KOSTEN[Varis_num,2])
  return(Total_costs)
}
