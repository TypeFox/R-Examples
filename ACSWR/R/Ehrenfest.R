Ehrenfest <-
function(n) {
  States <- c(0, seq(1,2*n))
  TPM <- matrix(0,nrow=length(States),ncol=length(States),dimnames=
                  list(seq(0,2*n),seq(0,2*n)))
  tran_prob <- function(i,n) {
    tranRow <- rep(0,2*n+1)
    if(i==0) tranRow[2] <- 1
    if(i==2*n) tranRow[(2*n+1)-1] <- 1
    if(i!=0 & i!=2*n) {
      j=i+1
      tranRow[j-1] <- i/(2*n)
      tranRow[j+1] <- 1-i/(2*n)
    }
    return(tranRow)
  }
  for(j in 0:(2*n))TPM[j+1,] <- tran_prob(j,n)
  return(TPM)
}
