EE <- function(MO,Output){
  #print("Entering EE")
  DesignLength=length(Output)
  VarO<-DesignLength*var(Output)/(DesignLength-1)
  I<-sum((Output-MO)^2)/DesignLength
  if (VarO==0 & I == 0) return(1)
  else return(1-I/VarO)
}