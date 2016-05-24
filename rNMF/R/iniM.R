
## Randomly initialize a matrix.

 ## Randomly initialize a matrix with uniform distribution.
 
 ## @param large The upper limit of the uniform distribution.
 ## @param nrow The number of rows of the matrix to be generated.
 ## @param ncol The number of columns of the matrix to be generated.
 ## @param small The lower limit of the uniform distribution.
 ## @param my.seed The seed to be used for the initialization of the matrix. 
 
 ## @return A nrow by ncol matrix.

 ## @examples
 ## initM(1,10,10)

initM = function(large, nrow, ncol, small = 0, my.seed = NULL){
    if(!missing(my.seed)) set.seed(my.seed)
    M = matrix(runif(nrow * ncol, small, large), nrow = nrow, ncol = ncol)
    return(M)
}
