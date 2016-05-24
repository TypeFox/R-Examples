descendants <- function(m,k){
  ## the done object indicates what rows of m were used
  done <- k
  if (m[k,1] < 0)
    left <- -m[k,1]
  else {
    junk <- descendants(m, m[k,1])
    left <- junk[[1]]
    done <- c(done, junk[[2]])
  }
  if (m[k,2] < 0)
    right <- -m[k,2]
  else {
    junk <- descendants(m, m[k,2])
    right <- junk[[1]]
    done <- c(done, junk[[2]])
  } 
  return(list(c(left, right), done))
}
