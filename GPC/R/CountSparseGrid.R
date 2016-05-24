CountSparseGrid <- function(d,L,Growth,Rule){
# 
# Input variables
#     d         : dumber of dimensions
#     L         : Quadrature level
#     Growth    :'Lin','LinOdd','LindEven','ExpOdd','ExpEven'
#     Rule      :'Gauss','Fejer','ClenshawCurtis','GaussKP' or 'LegendreKP','HermiteKP','HermiteKP126' or 'HermiteKP128'
  
  # rm(list=ls())
  # d = 3,L = 5,Growth = 'ExpOdd',Rule = 'ClenshawCurtis'
  
  MonomialSize    = GetNominalSize(L+1,Growth,Rule)#GetMonomialSize(L+1,Growth,Rule)
  StageSize       = rep(0,L)
  PCTable         = indexCardinal(d,L)
  
  for (P in 1:L){
      for (m in getM(d,P-2):(getM(d,P-1)-1)){StageSize[P] <- StageSize[P] + prod(MonomialSize$Unique[PCTable[,m+1]+1])}
      if (P != L){StageSize[P+1] <- StageSize[P]}
  }
  return(StageSize)
}