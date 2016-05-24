## compute sample entropy
##
## Input Parameters
##
## y  input signal vector
## M  maximum template length (default M=5)
## r  matching threshold (default r=.2)

## sflag  Standardize signal(default yes/sflag=1) 
## vflag  Calculate standard errors (default no/vflag=0) 


## Output Parameters
## e sample entropy estimates for m=0,1,...,M-1
## se standard error estimates for m=0,1,...,M-1
## A number of matches for m=1,...,M
## B number of matches for m=0,...,M-1

.packageName <- 'mousetrack'

sampen <- function(y, M, r, sflag, vflag){

## initialize the arguments with default values in case they are empty

    if (exists('M') == FALSE){ M =  5 }
    if (exists('r') == FALSE){ r = .2 }
    if (exists('sflag') == FALSE){ sflag = 1 }
    if (exists('vflag') == FALSE){ vflag = 0 }

    y = as.vector(y);
    n = length(y);

    if (sflag > 0){
        y = y - mean(y)
        s = sqrt(mean(y^2))   
        y = y/s
    
    }


    if (vflag > 0){
        se = sampense(y, M, r)
        se = se$se
    } else { se = vector()
         }


    res = as.data.frame( sampenc(y, M, r) )
    e = res$e

    return(list (se = se, e = e) )
}

# not sure why is the below repeated in matlab code

# N = n*(n-1)/2
# A = match(1:M)
# B = [N;A(1:(M-1))]
# N = n*(n-1)/2
