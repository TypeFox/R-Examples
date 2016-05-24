LOO <- function(MO,h,Output){
  #print("Entering LOO")
  DesignLength=length(Output)
  VarO<-DesignLength*var(Output)/(DesignLength-1)
  I<-sum(((Output-MO)/(1-h))^2)/DesignLength
  #print(h)
  #print(DesignLength)
  #print(I)
  #print(VarO)
  if (VarO==0 & I == 0) return(1)
  return(1-I/VarO)
}