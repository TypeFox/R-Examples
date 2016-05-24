
simulateGraph <- function(p,eta,extraeta=eta/5) {
  # ---------------------------------------------------------------
  # FUNCTION
  #   Draw at random a graph and create covariance matrix compatible
  #   with the graph.
  #   End-user function.
  # INPUT
  #   p: integer (number of nodes)
  #   eta: real number in (0,1) (proportion of edges in subgroups)
  #   extraeta: real number in (0,1) (proportion of edges itergroups)
  # OUTPUT
  #   MatG: p x p matrix (adjacency matrix of the graph)
  #   Nab: integer (number of nodes)
  #   mod: list, mod[[j]] = list of the neighbor of j
  #   C: p x p matrix (covariance matrix)
  #   PCor: p x p matrix (partial correlation matrix)
  # ---------------------------------------------------------------
  # Verify the arguments
   if (mode(p) != "numeric" || length(p) != 1 ||
       p <= 1 || (p%%1) !=   0) 
     stop("bad value of p")
   if (mode(eta) != "numeric" || length(eta) != 1 ||
       eta < 0 || eta > 1)
     stop("bad value of eta")
   if (mode(extraeta) != "numeric" || length(extraeta) != 1 ||
       extraeta < 0 || extraeta > 1)
     stop("bad value of extraeta")

   # Constants
   EPS <-  1/10
   PCORMIN  <-  1/1000
   PCORMAX  <- 5/1000
    N <- p*(p-1)/2
    Nab <- rbinom(1,N,extraeta)
    temp <- rep(0,N)
    nonzeros <- sample(1:N,Nab)
    temp[nonzeros] <- runif(length(nonzeros),min=-1,max=1)
    dimlab <- as.character(1:p)
    InfG <- matrix(0,nrow=p,ncol=p,
                   dimnames= list(dimlab, dimlab))
    InfG[lower.tri(InfG)] <- temp
#
    d <- floor(p/3)
    D <- d*(d-1)/2
    temp <- vector("integer", length=D)
    for (i in 1:3) {
    	Nab <- rbinom(1,D,eta)
   	temp[] <- 0
    	nonzeros <- sample(1:D,Nab)
    	temp[nonzeros] <- runif(length(nonzeros),min=-1,max=1)
    	G <- matrix(0,nrow=d,ncol=d)
    	G[lower.tri(G)] <- temp
	ind <- ((1+(i-1)*d):(i*d))
	InfG[ind,ind] <- G
    	}
#
    diag(InfG) <- runif(p,min=0,max=EPS)
    PCor <- InfG%*%t(InfG)/EPS
    PCor <- PCor+diag(runif(p,min= PCORMIN,max=PCORMAX))
    C <- solve(PCor)
    C <- cov2cor(C)
    PCor <- cov2cor(PCor)
    temp <- abs(sign(PCor))-diag(x=1,nrow=p,ncol=p)
    MatG <- matrix(0,nrow=p,ncol=p,
                   dimnames= list(dimlab, dimlab))
    MatG[lower.tri(MatG)] <- temp[lower.tri(temp)]
    Nab <- sum(MatG)

   # Build a list : each element j is the indexes of the elements
   # of temp[,j] that are equal to 1. It is zero, if none.
   mod <- apply(temp, 2, function(X, p) {
     r = (1:p)[X==1]
     if (length(r) == 0 ) r =0
     return(r)}, p)

    return(list(G=mattoSym(MatG),
                Nnodes=Nab,
                Neighb=modtoNeighb(mod),
                C=C,PCor=PCor))
} # fin simulateGraph

mattoSym <- function (G) {
  # ---------------------------------------------------------------
  # FUNCTION
  #   Transform G into a symmetric matrix, with diagonal equal to 0
  # INPUT
  #   G: matrix with nrow=ncol
  # OUTPUT
  #   G: idem as input, where the lower triangle is set equal
  #  to the upper one.
  # CALLED BY
  #   simulateGraph
  # ---------------------------------------------------------------
for (icol in 1:ncol(G)) {
  G[icol,] <- G[, icol]
  G[icol, icol] <- 0
}
return(G)
} # fin mattosym

modtoNeighb <- function (mod) {
  # ---------------------------------------------------------------
  # FUNCTION
  #   Transform the list of neighbours into a matrix
  # INPUT
  #   mod: list. Each component i is a vector with the
  #       neighbours of the vertex i. If none, scalar
  #       equal to zero
  # OUTPUT
  #   Neighb: matrix . The number of rows is the length of mod
  #       i.e, the number of vertices, and the number of columns
  #       is the maximal degree
  # CALLED BY
  #   simulateGraph
  # ---------------------------------------------------------------
  p <- length(mod)
  Neighb <- matrix(0, nrow=p, ncol=p)
  lmax <- 0
  for (imod in 1:p) {
    unmod <- mod[[imod]]
    if (unmod[1] != 0) {
      l <- 0
      for (a in unmod) {
        l <- l+1
        Neighb[imod,l] <- a
        if (l>lmax) lmax <- l
      } # fin a
    } # fin (unmod[1] != 0)
  } # fin imod
  Neighb <- matrix(Neighb[, 1:lmax], nrow=p)
  return(list(Neighb=Neighb, Dmax=lmax))
} # fin modtoNeighb
