forcerighthand<-function(U)
{
  ###  require(RFOC)
  ###  given an eigen vector decomposition,
  #   for the vectors to have a right handed system
  x12  = RSEIS::xprod(U[,1], U[,2])
  x3 = U[,3]
  if(!all(sign(x3)==sign(x12)))
    {
      U[,1] = -U[,1]
    }
  return(U)
}
