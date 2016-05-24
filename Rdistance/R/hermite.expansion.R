hermite.expansion <- function(x, expansions){
# Calculates hermite expansion for detection function.
# Input:
#       x = distances / w
#       expansions = number of expansion terms (1 - 5)
#
# Output:
#       expansion = a list of vectors, with 
#                   expansion[[1]] = expansion for the 1st term,
#                   expansion[[2]] = expansion for the 2nd, and so on.  
#
# Based on "probabilists' hermite polynomial as defined on Wikipedia.com

    if (expansions > 4){
        warning("Too many Hermite polynomial expansion terms. Only 4 used.")
        expansions = 4
    }
    
    if( expansions < 1 ) stop( "Number of expansions must be >= 1" )
  
   
    expansion = matrix(nrow=length(x), ncol=expansions)

    expansion[,1] = x^4 - 6*(x)^2 +3

    if(expansions >= 2){
        expansion[,2] = (x)^6 - 15*(x)^4 + 45*(x)^2 - 15
    } 
    if(expansions >= 3){
        expansion[,3] = (x)^8 - 28*(x)^6 + 210*(x)^4 - 420*(x)^2 + 105
    } 
    if(expansions >= 4){
        expansion[,4] =  (x)^10 - 45*(x)^8 + 630*(x)^6 - 3150*(x)^4 + 4725*(x)^2 - 945
    }         
    
    return(expansion)
            
}
