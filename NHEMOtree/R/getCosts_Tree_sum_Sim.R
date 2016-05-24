getCosts_Tree_sum_Sim <-
function(Kostenmatrix, Protein_IDs){
  Kostenvektor<- Kostenmatrix[Protein_IDs, 2]
  Gesamtkosten<- sum(Kostenvektor)
  return(Gesamtkosten)
}
