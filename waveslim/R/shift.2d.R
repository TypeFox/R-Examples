shift.2d <- function(z, inverse=FALSE) {
  ## "Center of Energy"
  coe <- function(g) { 
    sum(0:(length(g)-1) * g^2) / sum(g^2)
  }

  wf <- attributes(z)$wavelet
  h <- wave.filter(wf)$hpf
  g <- wave.filter(wf)$lpf
  
  J <- (length(z) - 1) / 3
  m <- nrow(z[[1]])
  n <- ncol(z[[1]])
  
  nu.H <- round(2^(1:J-1) * (coe(g) + coe(h)) - coe(g), 0)
  nu.Hm <- ifelse(nu.H/m < 1, nu.H, nu.H - trunc(nu.H/m) * m)
  nu.Hn <- ifelse(nu.H/n < 1, nu.H, nu.H - trunc(nu.H/n) * n)
  nu.G <- round((2^(1:J) - 1) * coe(g), 0)
  nu.Gm <- ifelse(nu.G/m < 1, nu.G, nu.G - trunc(nu.G/m) * m)
  nu.Gn <- ifelse(nu.G/n < 1, nu.G, nu.G - trunc(nu.G/n) * n)
  
  if (!inverse) {
    ## Apply the phase shifts
    for (j in 0:(J-1)) {
      Hm.order <- c((nu.H[j+1]+1):m, 1:nu.H[j+1])
      Hn.order <- c((nu.H[j+1]+1):n, 1:nu.H[j+1])
      Gm.order <- c((nu.G[j+1]+1):m, 1:nu.G[j+1])
      Gn.order <- c((nu.G[j+1]+1):n, 1:nu.G[j+1])
      z[[3*j+1]] <- z[[3*j+1]][Gm.order, Hn.order]
      z[[3*j+2]] <- z[[3*j+2]][Hm.order, Gn.order]
      z[[3*j+3]] <- z[[3*j+3]][Hm.order, Hn.order]
    } 
    z[[3*J+1]] <- z[[3*J+1]][Gm.order, Gn.order]
  } else {
    ## Apply the phase shifts "reversed"
    for (j in 0:(J-1)) {
      Hm.order <- c((m-nu.H[j+1]+1):m, 1:(m-nu.H[j+1]))
      Hn.order <- c((n-nu.H[j+1]+1):n, 1:(n-nu.H[j+1]))
      Gm.order <- c((m-nu.G[j+1]+1):m, 1:(m-nu.G[j+1]))
      Gn.order <- c((n-nu.G[j+1]+1):n, 1:(n-nu.G[j+1]))
      z[[3*j+1]] <- z[[3*j+1]][Gm.order, Hn.order]
      z[[3*j+2]] <- z[[3*j+2]][Hm.order, Gn.order]
      z[[3*j+3]] <- z[[3*j+3]][Hm.order, Hn.order]
    }
    z[[3*J+1]] <- z[[3*J+1]][Gm.order, Gn.order]
  }
  return(z)
}

