cplxdual2D <- function(x, J, Faf, af) {

  ## Dual-Tree Complex 2D Discrete Wavelet Transform
  ##
  ## USAGE:
  ##   w = cplxdual2D(x, J, Faf, af)
  ## INPUT:
  ##   x - 2-D array
  ##   J - number of stages
  ##   Faf{i}: first stage filters for tree i
  ##   af{i}:  filters for remaining stages on tree i
  ## OUTPUT:
  ##   w{j}{i}{d1}{d2} - wavelet coefficients
  ##       j = 1..J (scale)
  ##       i = 1 (real part); i = 2 (imag part)
  ##       d1 = 1,2; d2 = 1,2,3 (orientations)
  ##   w{J+1}{m}{n} - lowpass coefficients
  ##       d1 = 1,2; d2 = 1,2 
  ## EXAMPLE:
  ##   x = rand(256);
  ##   J = 5;
  ##   [Faf, Fsf] = FSfarras;
  ##   [af, sf] = dualfilt1;
  ##   w = cplxdual2D(x, J, Faf, af);
  ##   y = icplxdual2D(w, J, Fsf, sf);
  ##   err = x - y;
  ##   max(max(abs(err)))
  ##
  ## WAVELET SOFTWARE AT POLYTECHNIC UNIVERSITY, BROOKLYN, NY
  ## http://taco.poly.edu/WaveletSoftware/

  ## normalization
  x <- x/2
  w <- vector("list", J+1)
  
  for (m in 1:2) {
    w[[1]][[m]] <- vector("list", 2)
    for (n in 1:2) {
      w[[1]][[m]][[n]] <- vector("list", 2)
      temp <- afb2D(x, Faf[[m]], Faf[[n]])
      lo <- temp$lo
      w[[1]][[m]][[n]] <- temp$hi
      if (J > 1) {
        for (j in 2:J) {
          temp <- afb2D(lo, af[[m]], af[[n]])
          lo <- temp$lo
          w[[j]][[m]][[n]] <- temp$hi
        }
        w[[J+1]][[m]][[n]] <- lo
      }
    }
  }

  for (j in 1:J) {
    for (m in 1:3) {
      w[[j]][[1]][[1]][[m]] <- pm(w[[j]][[1]][[1]][[m]])
      w[[j]][[2]][[2]][[m]] <- pm(w[[j]][[2]][[2]][[m]])
      w[[j]][[1]][[2]][[m]] <- pm(w[[j]][[1]][[2]][[m]])
      w[[j]][[2]][[1]][[m]] <- pm(w[[j]][[2]][[1]][[m]])
    }
  }
  return(w)
}

icplxdual2D <- function(w, J, Fsf, sf) {

  ## Inverse Dual-Tree Complex 2D Discrete Wavelet Transform
  ## 
  ## USAGE:
  ##   y = icplxdual2D(w, J, Fsf, sf)
  ## INPUT:
  ##   w - wavelet coefficients
  ##   J - number of stages
  ##   Fsf - synthesis filters for final stage
  ##   sf - synthesis filters for preceeding stages
  ## OUTPUT:
  ##   y - output array
  ## See cplxdual2D
  ##
  ## WAVELET SOFTWARE AT POLYTECHNIC UNIVERSITY, BROOKLYN, NY
  ## http://taco.poly.edu/WaveletSoftware/

  for (j in 1:J) {
    for (m in 1:3) {
      w[[j]][[1]][[1]][[m]] <- pm(w[[j]][[1]][[1]][[m]])
      w[[j]][[2]][[2]][[m]] <- pm(w[[j]][[2]][[2]][[m]])
      w[[j]][[1]][[2]][[m]] <- pm(w[[j]][[1]][[2]][[m]])
      w[[j]][[2]][[1]][[m]] <- pm(w[[j]][[2]][[1]][[m]])
    }
  }

  y <- matrix(0, 2*nrow(w[[1]][[1]][[1]][[1]]), 2*ncol(w[[1]][[1]][[1]][[1]]))
  for (m in 1:2) {
    for (n in 1:2) {
      lo <- w[[J+1]][[m]][[n]]
      if (J > 1) {
        for (j in J:2) {
          lo <- sfb2D(lo, w[[j]][[m]][[n]], sf[[m]], sf[[n]])
        }
        lo <- sfb2D(lo, w[[1]][[m]][[n]], Fsf[[m]], Fsf[[n]])
        y <- y + lo
      }
    }
  }

  ## normalization
  return(y/2)
}

