##' transmission loss through a prism
##'
##' transmission loss through a prism
##' @title transmission
##' @export
##' @param n prism refractive index
##' @param external external incident angle in radians
##' @param polarisation polarisation
##' @return transmission 
##' @author baptiste Auguie
transmission <- function(n, external, polarisation = "p"){
  
  alpha <- asin(sin(external) / n) # % refracted angle
  
  if(polarisation == 'p'){
          4 * n*cos(external)*cos(alpha) /
            (n * cos(external) + cos(alpha))^2
        } else {
         .NotYetImplemented()
        }			
  
}
