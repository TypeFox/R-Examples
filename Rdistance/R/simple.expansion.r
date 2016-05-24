simple.expansion <- function(x, expansions){
# Calculates simple polynomial for detection function.
# Input:
#       x = distances / w
#       expansions = number of expansion terms (1 - 4)
#
# Output:
#       expansion = a list of vectors, with
#                   expansion[,1] = expansion for the 1st term,
#                   expansion[,2] = expansion for the 2nd, and so on.
#

    if (expansions > 4){
        warning("Too many Simple polynomial expansion terms. Only 4 used.")
        expansions = 4
    }

    if( expansions < 1 ) stop( "Number of expansions must be >= 1" )

    expansion = matrix(nrow=length(x), ncol=expansions)

    expansion[,1] = (x)^4

    if(expansions >= 2){
        expansion[,2] = (x)^6
    }
    if(expansions >= 3){
        expansion[,3] = (x)^8
    }
    if(expansions >= 4){
        expansion[,4] =  (x)^10
    }

    return(expansion)

}
