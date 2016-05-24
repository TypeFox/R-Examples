MOLS <- function(p, n, primpol = GF(p, n)[[2]][1, ]) {
  
  primen100 <- c(2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 
                 61, 67, 71, 73, 79, 83, 89, 97)
  if (!(p %in% primen100)) {
    stop("p is not a prime number less than 100.")
  }
  
  ## Now obsolete: require(conf.design) if( (p%%1) || (p < 2) ){stop('p is not a
  ## prime number.')} if( primes(p)[length(primes(p))]!=p ){stop('p is not a prime
  ## number.')}
  
  if (((n%%1) != 0) || (n < 1)) {
    stop("n is not a positive integer.")
  }
  
  ord <- p^n  # order of the field
  if ((p^n) < 3) {
    stop("The order of the field is too small")
  }
  
  M <- array(0, c(ord, ord, ord - 1))  # The resulting ord-1 latin squares are put into M
  
  if (n == 1) {
    ## ord is prime
    
    f <- 0:(p - 1)
    for (m in 1:(p - 1)) {
      for (i in 1:p) {
        for (j in 1:p) {
          dummy <- (f[m + 1] * f[i] + f[j])%%p  # a_ij(m) = f_m * f_i + f_j,
          ## i,j=0,1,..,ord-1; arithmetic is done modulo p
          M[i, j, m] <- dummy + 1
        }
      }
    }
  } else {
    ## ord is prime power
    f <- factor.comb(p, n)  # rows are polynomials f_0 through f_(ord-1)
    for (m in 1:(ord - 1)) {
      for (i in 1:ord) {
        for (j in 1:ord) {
          dummy <- mult(f[m + 1, ], f[i, ]) + c(numeric(n - 1), f[j, ])
          ## a_ij(m) = f_m * f_i + f_j, i,j=0,1,..,ord-1
          
          dummy <- (redu(dummy, primpol))%%p  # division mod p
          a <- 0
          for (r in 1:(n - 1)) {
            a <- a + dummy[r] * (p^(n - r))
          }
          M[i, j, m] <- a + dummy[n] + 1
          
        }
      }
    }
  }
  M
} 
