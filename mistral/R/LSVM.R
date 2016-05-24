#-----------------------------------------------------------------------------------------------------------------#
# Function LSVM
#-----------------------------------------------------------------------------------------------------------------#

#return an estimation of the class -1 or +1 of x
#' @export
LSVM  <- function(x, A.model.lsvm, convexity){

  if(is.null(dim(x))){
    rX  <- 1
    x   <- matrix(x, nrow = 1)
  }else{
    rX <- dim(x)[1]
  }
    
  DA <- dim(A.model.lsvm)
    
  if(is.null(DA)){
    rA  <- 1
    cA  <- rX + 1
    A <- matrix(x, nrow = 1)
  }else{ 
    rA <- DA[1]
    cA <- DA[2]
  }

  res <- apply(x, MARGIN = 1, function(u){  
                                temp <- numeric(rA)
                                for(j in 1:rA){
                                  temp[j] <- sign(A.model.lsvm[j, -cA]%*%u + A.model.lsvm[j, cA])
                                }
                                #test if u is below or above all the hyperplanes defined in the matrix A
                                res.temp <- ifelse( sum(temp) == convexity*rA, convexity, -convexity)
                                return(res.temp)
                              })    
  return(res)
}

