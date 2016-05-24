# Function to handle singulare matrices in GLARMA runs

mySolve <- function(A) {
  ErrCode <- 0
  Ainv <- try(solve(A), TRUE)
  if (any(is.na(Ainv))){
    ErrCode <- 1
    Ainv <- A
  }
  
  list(Ainv = Ainv, ErrCode = ErrCode)
} 
