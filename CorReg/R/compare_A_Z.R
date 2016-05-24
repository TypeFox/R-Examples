# ' description de la repartition des A vis a vis de la structure
compare_A_Z<-function (A = A, Z = Z) 
{
  if(length(A)==ncol(Z)+1){A=A[-1]}
  quiA = which(A != 0)
  quip2 = which(colSums(Z) != 0)
  p2 = length(quip2)
  nbA = length(quiA)
  nbAp2 = sum(duplicated(c(quiA, quip2)))
  quip1 = which(rowSums(Z) != 0)
  p1 = length(quip1)
  nbAp1 = sum(duplicated(c(quiA, quip1)))
  p3 = ncol(Z) - p1 - p2
  nbAp3 = length(which(A != 0)) - nbAp1 - nbAp2
  return(list(p2 = p2, p1 = p1, p3 = p3, nbAp2 = nbAp2, nbAp1 = nbAp1, 
              nbAp3 = nbAp3))
}