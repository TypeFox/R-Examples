dwpt.sim <- function(N, wf, delta, fG, M=2, adaptive=TRUE, epsilon=0.05) {
  M <- M*N
  J <- log(M, 2)
  jn <- rep(1:J, 2^(1:J))
  jl <- length(jn)

  if( adaptive ) {
    Basis <- find.adaptive.basis(wf, J, fG, epsilon) 
  } else {
    Basis <- numeric(jl)
    a <- min((1:jl)[jn == J])
    b <- max((1:jl)[jn == J])
    Basis[a:b] <- 1
  }

  Index <- (1:jl)[as.logical(Basis)]
  Length <- 2^jn

  variance <- bandpass.var.spp(delta, fG, J, Basis, Length)
  z <- vector("list", jl)
  class(z) <- "dwpt"
  attr(z, "wavelet") <- wf

  for(i in Index)
    z[[i]] <- rnorm(M/Length[i], sd=sqrt(Length[i]*variance[i]))

  x <- idwpt(z, Basis)
  xi <- trunc(runif(1, 1, M-N))
  return(x[xi:(xi+N-1)])
}

find.adaptive.basis <- function(wf, J, fG, eps) {
  H <- function(f, L) {
    H <- 0
    for(l in 0:(L/2-1))
      H <- H + choose(L/2+l-1,l) * cos(pi*f)^(2*l)
    H <- 2 * sin(pi*f)^L * H
    return(H)
  }

  G <- function(f, L) {
    G <- 0
    for(l in 0:(L/2-1))
      G <- G + choose(L/2+l-1,l) * sin(pi*f)^(2*l)
    G <- 2 * cos(pi*f)^L * G
    return(G)
  }
  
  L <- wave.filter(wf)$length
  jn <- rep(1:J, 2^(1:J))
  jl <- length(jn)
  U <- numeric(jl)
  U[1] <- G(fG, L)
  U[2] <- H(fG, L)
  for(j in 2:J) {
    jj <- min((1:jl)[jn == j])
    jp <- (1:jl)[jn == j-1]
    for(n in 0:(2^j/2-1)) {
      if (n%%2 == 0) {
        U[jj + 2 * n + 1] <- U[jp[n+1]] * H(2^(j-1)*fG, L)
        U[jj + 2 * n] <- U[jp[n+1]] * G(2^(j-1)*fG, L)
      } else {
        U[jj + 2 * n] <- U[jp[n+1]] * H(2^(j-1)*fG, L)
        U[jj + 2 * n + 1] <- U[jp[n+1]] * G(2^(j-1)*fG, L)
      }
    }
  }
  return(ortho.basis(U < eps))
}

bandpass.var.spp <- function(delta, fG, J, Basis, Length) {
  a <- unlist(sapply(2^(1:J)-1, seq, from=0, by=1)) / (2*Length)
  b <- unlist(sapply(2^(1:J), seq, from=1, by=1)) / (2*Length)
  bp.var <- rep(0, length(Basis))
  for(jn in (1:length(Basis))[as.logical(Basis)]) {
    if(fG < a[jn] | fG > b[jn])
      bp.var[jn] <- 2*integrate(spp.sdf, a[jn], b[jn], d=delta, fG=fG)$value
    else {
      result1 <- 2*integrate(spp.sdf, a[jn], fG, d=delta, fG=fG)$value
      result2 <- 2*integrate(spp.sdf, fG, b[jn], d=delta, fG=fG)$value
      bp.var[jn] <- result1 + result2
    }
  }
  return(bp.var)
}

