generatePCEcoeff <- function(m,n,y,PhiZn,weights){
  PCEcoeff = rep(0,m)
  PhiIJ = rep(0,m)
  # print(PhiZn)
  for (mm in 1:m){
    for (zz in 1:n){
      PCEcoeff[mm] = PCEcoeff[mm] + y[zz]*PhiZn[mm,zz]*weights[zz] 
      PhiIJ[mm] = PhiIJ[mm] + PhiZn[mm,zz]*PhiZn[mm,zz]*weights[zz]
    }
    PCEcoeff[mm] = PCEcoeff[mm]/PhiIJ[mm]
  }
  res = list(PCEcoeff = PCEcoeff, PhiIJ = PhiIJ)
  return(res)
}
