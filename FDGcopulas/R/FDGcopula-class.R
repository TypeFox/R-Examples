
## FDGcopula object 
setClass(
  Class="FDGcopula",
  representation=representation(
    dimension = "integer",
    parameters = "numeric",
    family = "character",
    extremevalue = "logical",
    parameterrange = "numeric")
  )
  


  
## Constructor 
FDGcopula <- function(family, parameters, extremevalue=FALSE, checkbounds=TRUE){
  copula <- new(Class="FDGcopula")
  copula@parameters <- parameters
  copula@dimension <- length(parameters)
  copula@family <- family
  copula@extremevalue <- extremevalue
  parameterRange <- switch(copula@family,
                           "frechet" = c(0,1),
                           "cuadrasauge" = c(0,1),
                           "sinus" = c(0,pi/2),
                           "exponential" = c(0,Inf))
  theta <- copula@parameters
  copula@parameterrange <- parameterRange
  if(checkbounds){
      if(sum(theta>parameterRange[2] | theta<parameterRange[1])>=1){
          stop("The parameters lie outside their admissible space.")
      }else{}
  }else{
      if(length(copula@parameters)!=copula@dimension){
          stop("Dimension should be equal to the number of parameters")
          }else{}
   }
  copula
}
 
 




