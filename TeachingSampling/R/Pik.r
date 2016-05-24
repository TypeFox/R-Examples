Pik <- function(p, Ind){
multip <- p*Ind
pik <- colSums(multip)
t(pik)
}