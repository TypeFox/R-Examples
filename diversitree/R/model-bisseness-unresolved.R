## Compute the probabilities of generating every state, given that we
## started with one species in state 0 or state 1.  For lt = length(t)
## times, returns a n(space) x lt x 2 array.
##   'scal': a multiplier to get use of extra precision (the vector
##           will be scaled so that the total probability sums to
##           scal.  Underflow is possible where scal=1, though).
##   'tol': tolerance, used as exit condition (when changes in the
##          probability space do not exceed tol, the step size is good
##          enough).  Lower values _may_ be faster and _may_ have
##          lower accuracy (but with large numbers of clades being
##          calculated, this effect can be small).  It is not clear to
##          me if tol should be resacled by scal before being passed
##          in.
##   'm': parameter affecting internal arguments.  Small numbers often
##        faster, sometimes slower.  Play around.
nucexp <- function(nt, la0, la1, mu0, mu1, q01, q10, p0c, p0a,
                   p1c, p1a, t, scal=1, tol=1e-10, m=15) {
  stopifnot(all(c(la0, la1, mu0, mu1, q01, q10, p0c, p0a, p1c, p1a) >= 0))
  lt <- length(t)
  stopifnot(lt > 0)
  stopifnot(scal >= 1)
  stopifnot(m > 1)
  n <- (nt*(nt+1)/2+1)
  res <- .Fortran("nucexp",
                  nt   = as.integer(nt),
                  la0  = as.numeric(la0),
                  la1  = as.numeric(la1),
                  mu0  = as.numeric(mu0),
                  mu1  = as.numeric(mu1),
                  q01  = as.numeric(q01),
                  q10  = as.numeric(q10),
                  p0c  = as.numeric(p0c),
                  p0a  = as.numeric(p0a),
                  p1c  = as.numeric(p1c),
                  p1a  = as.numeric(p1a),
                  t    = as.numeric(t),
                  lt   = as.integer(lt),
                  scal = as.numeric(scal),
                  tol  = as.numeric(tol),
                  m    = as.integer(m),
                  w    = numeric(2*n*lt),
                  flag = integer(1))
  if ( res$flag < 0 )
    stop("Failure in the Fortran code")
  else if ( res$flag > 0 )
    cat(sprintf("Flag: %d\n", res$flag))
  array(res$w, c(n, lt, 2))
}

## This is more useful; given model parameters, and vectors of times
## (t), number of species in each clade (Nc), numbers of species in
## known states (nsc) and numbers of species in state b (k), calculate
## the probabilities of generating the clade, and of extinction.
## Returns a length(t) x 4 matrix, where the first two columns are the
## probabilities of generating the clade, and the second two are the
## probabilities of extinction.
nucexpl <- function(nt, la0, la1, mu0, mu1, q01, q10, p0c, p0a,
                    p1c, p1a, t, Nc, nsc, k, scal=1, tol=1e-10, m=15) {
  tt <- sort(unique(t))
  ti <- as.integer(factor(t))

  lt <- length(tt)
  lc <- length(Nc)
  stopifnot(lc > 0, lt > 0, max(Nc) < nt,
            length(Nc) == lc, length(nsc) == lc, length(k) == lc,
            all(nsc <= Nc), all(nsc >= 0), all(k <= nsc),
            scal >= 1, m > 1)

  res <- .Fortran("nucexpl",
                  nt   = as.integer(nt),
                  la0  = as.numeric(la0),
                  la1  = as.numeric(la1),
                  mu0  = as.numeric(mu0),
                  mu1  = as.numeric(mu1),
                  q01  = as.numeric(q01),
                  q10  = as.numeric(q10),
                  p0c  = as.numeric(p0c),
                  p0a  = as.numeric(p0a),
                  p1c  = as.numeric(p1c),
                  p1a  = as.numeric(p1a),
                  t    = as.numeric(tt),
                  lt   = as.integer(lt),
                  ti   = as.integer(ti),
                  Nc   = as.integer(Nc),
                  nsc  = as.integer(nsc),
                  k    = as.integer(k),
                  lc   = as.integer(lc),
                  scal = as.numeric(scal),
                  tol  = as.numeric(tol),
                  m   = as.integer(m),
                  ans  = numeric(4*lc),
                  flag = integer(1))
  if ( res$flag < 0 )
    stop("Failure in the Fortran code")
  else if ( res$flag > 0 )
    cat(sprintf("Flag: %d\n", res$flag))
  matrix(res$ans, ncol=4)
}

## Construct the transition matrix.  Again, non-R style, as this was
## used as a template for constructing the same in Fortran.
make.matrix.ness <- function(nt, mu.a, mu.b, lambda.a, lambda.b, q.ba,
                        q.ab, p.ac, p.aa, p.bc, p.ba) {
  ## Diagonals (are linear sums of these: -(na*pa + nb*pb))
  ## The chance of no change is unaffected by BiSSEness, which only alters what
  ## happens at speciation.

  pa <- lambda.a + mu.a + q.ba
  pb <- lambda.b + mu.b + q.ab

  ## A few useful numbers:
  ## n1 represents the state space size, excluding the absorbing state and the state with nt-1 species
  ## 2*n1 also represents the total number of non-zero elements in each of the transition and extinction matrices
  ## nm represents the state space size:  {0,0}, {1,0}, {0,1}, ...{absorbing with nt species}
  n1 <- nt*(nt-1)/2
  nm <- n1 + nt + 1

  ## Create some useful vectors
  ## j1 contains the positions in the vector space that have some species in state a (na>0)
  ## j1+1 contains the positions in the vector space that have some species in state b (nb>0)
  ## na contains just the vector of numbers of species in state a (dropping those with na=0)
  ## nb contains just the vector of numbers of species in state b (dropping those with nb=0)
  j1 <- integer(n1)
  na <- integer(n1)
  nb <- integer(n1)
  k <- 1
  for ( i in 1:(nt-1) ) {
    for ( j in 1:i ) {
      j1[k] <- i + k
      na[k] <- i - j + 1
      nb[k] <- j
      k <- k + 1
    }
  }

  ## Character state transitions and extinctions
  m <- matrix(0, nm, nm)
  for ( i in 1:n1 ) {
    m[j1[i],   i]       <- na[i]*mu.a
    m[j1[i]+1, i]       <- nb[i]*mu.b
    m[j1[i],   j1[i]+1] <- na[i]*q.ba
    m[j1[i]+1, j1[i]]   <- nb[i]*q.ab
  }

  ## Speciation and diagonals (BiSSE):
##  for ( i in 2:n1 ) {
##    if ( na[i] > 1 )
##      m[i, j1[i]]   <- (na[i]-1)*lambda.a
##    if ( nb[i] > 1 )
##    m[i, j1[i]+1] <- (nb[i]-1)*lambda.b
##    m[i,i]        <- - (na[i]-1)*pa - (nb[i]-1)*pb
##  }

  ## Speciation and diagonals (BiSSEness):
  ## (Calculations are performed where na[i]-1 represents number of state a species in source.)
  for ( i in 2:n1 ) {
    if ( na[i] > 1 ) {
      m[i, j1[i]]   <- (na[i]-1)*lambda.a*(1-p.ac)
      m[i, j1[i]+1]   <- (na[i]-1)*lambda.a*p.ac*p.aa
      m[i, j1[i]+2]   <- (na[i]-1)*lambda.a*p.ac*(1-p.aa)
    }
    if ( nb[i] > 1 ) {
      m[i, j1[i]+1] <- m[i, j1[i]+1] + (nb[i]-1)*lambda.b*(1-p.bc)
      m[i, j1[i]] <- m[i, j1[i]] + (nb[i]-1)*lambda.b*p.bc*p.ba
##    While the above may be contributed by speciation in state a, the following cannot be.  
      m[i, j1[i]-1] <- (nb[i]-1)*lambda.b*p.bc*(1-p.ba)
    }
    m[i,i]        <- - (na[i]-1)*pa - (nb[i]-1)*pb
  }

  ## Speciation in the special final column, diagonals for the last
  ## class (unaffected by bisseness, because all movements to last
  ## class are collapsed, regardless of trait combination).
  k <- nt*(nt-1)/2
  for ( i in 1:nt ) {
    m[k+i, nm]  <-  (nt-i)*lambda.a + (i-1)*lambda.b
    m[k+i, k+i] <- -(nt-i)*pa       - (i-1)*pb
  }

  m
}
