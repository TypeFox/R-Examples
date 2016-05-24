fillMatrix <-
function(totalMatrix, smallerMatrix) {
  ## loop for the smaller matrix
  for(k in 1:length(smallerMatrix[1,])) {
    for(l in 1:length(smallerMatrix[,1])) {
      counter = as.integer(0)
      ## i and j are indices for the bigger matrix
      ## counting elementwise for 3x3 matrices the number of zeros
      i = 3*k-2;j = 3*l-2
      ## limits for the for loops of the small 3x3 matrices
      a = i+2; b = j+2
      if(i <= length(totalMatrix[1,])) {
        for(i in i:a) {
          for(j in j:b){
            if(totalMatrix[i,j] == 1) {
              counter = counter + 1
            } # end if
            ## reset j 
            j = 3*l-2
          } # end for j
        } # end for i
        ## if at least one of the elements of the matrix is 1, 
        ## then value for smaller matrix is 1 otherwise 0
        if(counter > 0) {
          smallerMatrix[k,l] = 1
          counter = as.integer(0)
        } # end if
        else {
          smallerMatrix[k,l] = 0
          counter = as.integer(0)
        } # end else
      } # end if
    } # end for l
  } # end for k
  return(smallerMatrix)
}
