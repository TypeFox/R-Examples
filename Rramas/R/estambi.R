estambi <-function (mat, matsd, equalsign) 
{
   mat<- as.matrix(mat)
   matsd <- as.matrix(matsd)
    if (equalsign != TRUE) {
        mat <- matrix(
                      rnorm(length(mat), mean = as.vector(mat), sd = as.vector(matsd)), 
            nrow = dim(mat)[1], ncol = dim(mat)[2])
    }
    if (equalsign == TRUE) {
      #get normal deviates (mean 0, sd = matsd), all with the same random sign
       deviates<- abs(rnorm(length(mat), mean = 0, sd = as.vector(matsd))) * sample(c(-1,1),1)
        mat <- mat + matrix(deviates, nrow = dim(mat)[1], ncol = dim(mat)[2])
    }
    # correct negative transitions and probabilities > 1
    mat[mat < 0] <- 0
    mat[-1,][mat[-1,] > 1] <- 1
    return(mat)
}

