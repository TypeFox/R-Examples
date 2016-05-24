CDF <- function(c, density){
     resultCFD <- integrate(density, lower = c, upper = Inf)$value
     return(resultCFD)
}









#