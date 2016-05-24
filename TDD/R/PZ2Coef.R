PZ2Coef = function(PZ, dt){
  # takes a poles/zeros object and returns b, a coefficients for filtering.
  # note that this is just an approximation for nonzero dt!
#  require(pracma)
  if(length(PZ$Knorm) == 0 ){
    PZ$Knorm = 1
  }
  if(length(PZ$Sense) == 0){
    PZ$Sense = 1
  }
  # calculate numerator polynomial (b is moving-average coefficients)
  b = ConvDiffCoef(rev(Poly(PZ$zeros)), dt)
  
  # calculate denominator polynomial (a is auto-regressive coefficients)
  a = ConvDiffCoef(rev(Poly(PZ$poles))/(PZ$Knorm * PZ$Sense), dt)

  return(list(b = b, a = a))
}

