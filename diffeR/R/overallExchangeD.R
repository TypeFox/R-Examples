overallExchangeD <- function(ctmatrix){
Exchange <- exchangeDj(ctmatrix)
overallexcd <- sum(Exchange)/2
return(overallexcd)
}
