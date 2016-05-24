# generate all possible haplotypes for n biallelic markers
genhaplopairs <- function(n) {
  two2n <- 2^n
# initiate and fill haplotype matrix bmat
  bmat <- matrix(0,two2n,n)
  bmat[1,1] <- 0
  bmat[2,1] <- 1
  two2i <- 1
  for (i in 2:n) {
    two2i <- two2i*2
    bmat[1:two2i,i] <- 0
    bmat[two2i+(1:two2i),i] <- 1
    bmat[two2i+(1:two2i),1:(i-1)] <- bmat[1:two2i,1:(i-1)]
  }
  np <- sum(choose(n,1:n)*2^(0:(n-1)))
  haplopairs <- matrix(0, np, 2)
  nstart <- 1
  nend <- 0
  g1indextable <- matrix(0,two2n,2)
  for (i in 2:two2n) {
    g1indextable[i,1] <- nstart
    j <- sum(bmat[i,])
    two2j <- 2^(j-1)
    g1indextable[i,2] <- two2j
    nend <- nend + two2j
    jj <- which(bmat[i,]==1) - 1
    haplopairs[nstart:nend, 1] <- bmat[1:two2j,1:j]%*%2^jj
    haplopairs[nstart:nend, 2] <- sum(2^jj) - haplopairs[nstart:nend, 1]
    nstart <- nstart + two2j
  }
  list(g1tbl=g1indextable,hpair=haplopairs)
}

