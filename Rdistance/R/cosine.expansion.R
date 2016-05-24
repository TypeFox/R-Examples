cosine.expansion <- function(x, expansions){
# Calculates cosine expansion for detection function.
# Input:
#       x = distances / w
#       expansions = number of expansion terms (1 - 5)
#
# Output:
#       expansion = a matrix with columns
#                   expansion[,1] = expansion for the 1st term,
#                   expansion[,2] = expansion for the 2nd, and so on.  
# changed 1st coeff to 2 - even function - jg
   
    if (expansions > 5){
        warning("Too many Cosine polynomial expansion terms. Only 5 used.")
        expansions = 5
    }
    
    if( expansions < 1 ) stop( "Number of expansions must be >= 1" )
    
    expansion = matrix(nrow=length(x), ncol=expansions)

    expansion[,1] = cos(2*pi*x)
    
    # I realize I could do this in a for loop, but I think this is faster.    
    if(expansions >= 2){
        expansion[,2] = cos(3*pi*x)
    }
    if(expansions >= 3){
        expansion[,3] = cos(4*pi*x)
    }
    if(expansions >= 4){
        expansion[,4]=  cos(5*pi*x)
    }
    if(expansions >= 5){
        expansion[,5] = cos(6*pi*x)
    }      
   
    return(expansion)
            
}
