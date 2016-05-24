LOOLar <- function(MO,h,Output,InvIM,dim,p){
  #print("Entering LOO")
  DesignLength=length(Output)
  VarO<-DesignLength*var(Output)/(DesignLength-1)
  I<-sum(((Output-MO)/(1-h))^2)/DesignLength
  P=choose(dim+p,p)
  CorrectFactor=1
  if(P<DesignLength){CorrectFactor=DesignLength*(1+sum(diag(InvIM)))/(DesignLength-P)}
  #print(h)
  #print(DesignLength)
  #print(I)
  #print(VarO)
  if (VarO==0 & I == 0) return(1)
  return(1-CorrectFactor*I/VarO)
}
