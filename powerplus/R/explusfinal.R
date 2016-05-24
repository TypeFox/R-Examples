explus <- function(a, numer, denom = 1){
  # Raises any base (real or complex) to any real power.
  #
  # Args:
  #   a: any base (real or complex)
  #   numer: numerator of (rational) exponent. Can be a decimal.
  #   denom: denominator of rational exponent (1 by default)
  #
  # Returns:
  #   The solution to the exponentiation operation supplied
  #   Returns a real-valued root whenever possible
  #   Returns the (principal) complex root if that is the only root
  if (numer %% 1 == 0){
    vec <- reduce.fraction (c(numer, denom))
    numer <- vec[1]
    denom <- vec[2]
    if (is.complex(a) == T){
      a ^ (numer / denom)
    }else{
      if (a < 0){
        if (numer %% 2 == 0){
          (a ^ numer) ^ (1 / denom)
        }else{
          if (denom %% 2 == 0){
            (as.complex(a)) ^ (numer / denom)
          }else{
            Nthroot <- function(x, n){ #only for n odd
              sign(x) * abs(x) ^ (1 / n)
            }
            (Nthroot(a, denom)) ^ numer
          }
        }
      }else{
        a ^ (numer / denom)
      }
    } 
  }else{
    denom <- as.numeric(unlist(strsplit(attributes(fractions(numer))$fracs,split="/")))[2]
    numer <- as.numeric(unlist(strsplit(attributes(fractions(numer))$fracs,split="/")))[1]
    if (is.complex(a) == T){
      a ^ (numer / denom)
    }else{
      if (a < 0){
        if (numer %% 2 == 0){
          (a ^ numer) ^ (1 / denom)
        }else{
          if (denom %% 2 == 0){
            (as.complex(a)) ^ (numer / denom)
          }else{
            Nthroot <- function(x, n){ #only for n odd
              sign(x) * abs(x) ^ (1 / n)
            }
            (Nthroot(a, denom)) ^ numer
          }
        }
      }else{
        a ^ (numer / denom)
      }
    }
  }    
}