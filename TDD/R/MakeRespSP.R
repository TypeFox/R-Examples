MakeRespSP = function(h, omega0, Sense, f_Sense = NULL){
  PZ = list(Sense = Sense,
    poles = -h * omega0 + c(-1, 1) * omega0 * sqrt(h^2 - 1+0i),
    np = 2,
    zeros = c(0, 0),
    nz = 2)
  
  if(is.null(f_Sense)){
    PZ$Knorm = 1
  }else{
    PZ$Knorm = 1/abs((2i*pi*f_Sense)^2/prod(2i*pi*f_Sense - PZ$poles))
  }
  PZ
}
  
