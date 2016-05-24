Matpow <- function(M, numer, denom = 1){
  # Raises any diagonalizable Matrix to any real power
  # Eigenvectors can be complex
  # Capabilities for non-diagonalizable matrices:
  #   Exponentiation (integer powers only)
  #   Square root
  #
  # Args:
  #   M: a square Matrix
  #   numer: numerator of (rational) exponent. Can be a decimal.
  #   denom: denominator of rational exponent (1 by default)
  #
  # Returns:
  #   The solution to the exponentiation operation supplied
  #   Allows to take powers and roots
  #   Allows to compute the inverse matrix (if invertible)
  #   Method based on spectral decomposition
  #   Returns a real-valued root whenever possible
  #   Returns the (principal) complex root if that is the only root
  frac <- (numer / denom)
  e.val <- eigen(M)$values
  n <- length(e.val)
  V <- t(eigen(M)$vectors)
  D <- diag(e.val)
  dexp <- numeric(n)
  if (rankMatrix(V) != n){
    if(frac == 0.5){
      sqrtm(M)
    }else{
      if(numer %% 1 == 0){
        M %^% numer
      }else{
        print("Sorry, requested operation requires a diagonalizable matrix.") 
      }
    }
  }else{
    if (numer %% 1 == 0){
      vec <- reduce.fraction (c(numer, denom))
      numer <- vec[1]
      denom <- vec[2]
      for (i in 1:n){  #allows for complex entries
        if (is.complex(e.val[i]) == T){
          dexp[i] <- e.val[i] ^ (numer / denom)
        }else{
          if (e.val[i] < 0){
            if (numer %% 2 == 0){
              dexp[i] <- ((e.val[i] ^ numer) ^ (1 / denom))
            }else{
              if (denom %% 2 == 0){
                dexp[i] <- ((as.complex(e.val[i])) ^ (numer / denom))
              }else{
                Nthroot <- function(x, n){ #only for n odd
                  sign(x) * abs(x) ^ (1 / n)
                }
                dexp[i] <- ((Nthroot(e.val[i], denom)) ^ numer)
              }
            }
          }else{
            dexp[i] <- (e.val[i] ^ (numer / denom))
          }
        }
      }  
    }else{
      denom <- as.numeric(unlist(strsplit(attributes(fractions(numer))$fracs,split="/")))[2]
      numer <- as.numeric(unlist(strsplit(attributes(fractions(numer))$fracs,split="/")))[1]
      for (i in 1:n){
        if (is.complex(e.val[i]) == T){
          dexp[i] <- e.val[i] ^ (numer / denom)
        }else{
          if (e.val[i] < 0){
            if (numer %% 2 == 0){
              dexp[i] <- ((e.val[i] ^ numer) ^ (1 / denom))
            }else{
              if (denom %% 2 == 0){
                dexp[i] <- ((as.complex(e.val[i])) ^ (numer / denom))
              }else{
                Nthroot <- function(x, n){
                  sign(x) * abs(x) ^ (1 / n)
                }
                dexp[i] <- ((Nthroot(e.val[i], denom)) ^ numer)
              }
            }
          }else{
            dexp[i] <- (e.val[i] ^ (numer / denom))
          }
        }  
      }
    }
    t(solve(V, dexp * V))
  }
}