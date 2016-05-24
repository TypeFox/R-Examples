# function used within calc.r to calculate locus-wise Hs values
Hs.value.calculation <- function(Hj.values2.1table){
Hs.one.locus <- Hs(as.numeric(as.vector(Hj.values2.1table$Hj.value)))
Hs.one.locus
}
