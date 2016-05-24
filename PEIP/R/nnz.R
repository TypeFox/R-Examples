nnz <-
function(h)
  {
    ###   return number of nonb-zero elements in a vector
    ###  duplicates a matlab function
    return(length(which(h>0)))
  }
