Sim_Costs <-
function(){
  Names <- c("X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8", "X9","X10")
  Kosten<- c(5,5,5,4,3,2,1,1,1,1)
  
  # Scaling costs
  Kosten<- Kosten/sum(Kosten)*100
  
  # Output
  KOSTEN<- as.data.frame(cbind(Names, as.data.frame(Kosten)))
  return(KOSTEN)
}
