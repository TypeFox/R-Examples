## Runif <- function(nrEdges, lB=0.1, uB=1) {
##   runif(nrEdges, lB, uB)
## }
m2g <- function(m) {
    ## INPUT: valid adjacency matrix
    ## m[i, j] = 1 corresponds to i -> j with weight 1
    ## OUTPUT: corresponding graphNEL object
    ## Improves 'as(Q, "graphNEL")' since neg weights are now OK
    ## m2g is tested in bugs: Jun15_randDAG
    n <- ncol(m)
    V <- colnames(m)
    edL <- vector("list", n)
    nmbEdges <- 0L
    for (i in seq_len(n)) {
        idx <- which( m[i, ] != 0 )
        nmbEdges <- nmbEdges + length(idx)
        edL[[i]] <- list(edges = idx, weights = m[i, idx])
    }
    if (nmbEdges > 0) {
	names(edL) <- V
	new("graphNEL", nodes = V, edgeL = edL, edgemode = "directed")
    }
    else
	new("graphNEL", nodes = V, edgemode = "directed")
}

## construct DAG (matrix Q) recursively ###########################
sampleQ <- function(n, K, p.w=1/2) {

  Q <- matrix(0, n, n)
  I <- length(K)
  j <- K[I]
  if(I >= 2) for(i in I:2) { # i = I, I-1, ..., 2
    Ki_1 <- K[i-1L]
    j. <- j + Ki_1
    for(p in (j-K[i]+1):j) {
      repeat{
        for(m in (j+1):j.) { ## FIXME rbinom(1, nn, *)
          Q[m, p] <- rbinom(1, 1, p.w)
        }
        if(sum(Q[(j+1):j., p])>0)
          break
      }
      if(i > 2) {
        for(m in (j.+1):n) { ## FIXME rbinom(1, nn, *)
          Q[m, p] <- rbinom(1, 1, p.w)
        }
      }
    }
    j <- j.
  }
  Q
}

################################################## ---> ../man/randDAG.Rd
##							-----------------
##' Random DAGs from igraph
randDAG <- function(n, d, method = "er", par1=NULL, par2=NULL, DAG = TRUE,
                    weighted = TRUE, wFUN = list(runif, min=0.1, max=1))
{
  if(!is.list(wFUN))
    wFUN <- list(wFUN)

  ## From igraph graph to upper triangular Q - via random re-labeling
  g2Q <- function(g, sparse=FALSE) {
    Q <- get.adjacency(g, sparse=sparse)
    perm <- sample.int(n)
    Q <- Q[perm, perm]
    Q * upper.tri(Q)
  }

  switch(method,
         "er" = {
           Q <- erDAG(n, p = d / (n-1))
         },

         "regular" = {
           ## s = d = number of neighbours in expectation
           s <- d
           g <- k.regular.game(n, s)
	   Q <- g2Q(g)
         },

         "watts" = {
           ## par1 = beta = fraction of interpolating between regular lattice (0) and ER (1)
           beta <- if(is.numeric(par1)) par1 else 1/2
           ## s = number of neighbours in expectation:
           s <- round(d/2)
           g <- watts.strogatz.game(1, n, s, beta)
	   Q <- g2Q(simplify(g))
         },

         "bipartite" = {
           ## par1 = alpha := fraction of part one
           alpha <- if(is.numeric(par1)) par1 else 1/2
           ## p = probability of connecting between two parts :
           p <- d/(2*alpha*(1-alpha)*n)
           n1 <- ceiling(n*alpha)
           n2 <- floor(n*(1-alpha))
           g <- bipartite.random.game(n1=n1, n2=n2, type="gnp", p=p)
	   Q <- g2Q(g)
         },

         "barabasi" = {
           ## par1 = power = power of preferential attachment
           power <- if(is.numeric(par1)) par1 else 1
           ## m = number of nodes adding in each step
           m <- round(d/2)

           mbar <- (2*m*n-m^2+m)/(2*(n-m))
           seq <- sample(c(m, m+1), n, replace=TRUE, prob=c(1-mbar+m, mbar-m))

           g <- barabasi.game(n=n, power=power, out.seq=seq, out.pref=TRUE, directed=FALSE)
	   Q <- g2Q(simplify(g))
         },

         "geometric" = {
           ## par1: dim=dimension:
           ## par2: if "geo", then weights are reciprocal to distances
           dim <- if(is.numeric(par1)) par1 else 2
           geo2 <- (is.character(par2) &&  par2 == "geo") # T / F
           ## r = euclidian radius :
           r <- ( (d*gamma(dim/2+1)) / ((n-1)*pi^(dim/2)) )^(1/dim)
           return( geoDAG(n, r, dim=dim, geo=geo2, DAG=DAG, weighted=weighted, wFUN=wFUN) )
         },

         "power" = {
           gamma <- findGamma(n, d)
           g <- powerLawDAG(n, gamma)
           Q <- g2Q(g)
         },

         "interEr" = {
           ## par1 = s.island = number of islands, n/s.island should be a positive integer
           ## par2 = alpha = fraction of inter connectivity

           s.island <- if(is.numeric(par1)) par1 else 2
           alpha    <- if(is.numeric(par2)) par2 else 1/4
           stopifnot(n %% s.island == 0)

           ## p = probability of connecting edges intra
           p <- 2*d*s.island / ((n-s.island)*(2+alpha))
           stopifnot(0 < p, p <= 1)

           m <- round((alpha*p*(n-2*s.island)*(n-s.island))/(2*s.island^2*(s.island-1)))
           g <- interconnected.islands.game(s.island, n/s.island, p, m)
	   Q <- g2Q(g)
         },
         stop("unsupported 'method': ", method))## switch end

  ## return
  undirunweight.to.dirweight(Q, n, DAG=DAG, weighted=weighted, wFUN=wFUN)
}

## AUX-FUNCTIONS ####################################
undirunweight.to.dirweight <- function(Q, n, DAG, weighted,
                                       wFUN = list(runif, min=0.1, max = 1)) {
## input Q: upper triangular matrix
  if(weighted) {
    nrEdge <- sum(Q)
    Q[Q==1] <- do.call(wFUN[[1]], c(nrEdge, wFUN[-1]))
  }

  if(!DAG) Q <- Q+t(Q)

  perm <- sample.int(n)
  Q <- Q[perm, perm]
  colnames(Q) <- rownames(Q) <- 1:n
  ## as(Q, "graphNEL")
  m2g(Q)
}

powerLawDAG <- function(n, gamma, maxtry = 20L) {
  ## generate power-law distribution
  stopifnot((n <- as.integer(n)) >= 2,
            (maxtry <- as.integer(maxtry)) >= 1)
  in1 <- seq_len(n - 1L)
  dist <- in1^(-gamma)
  dist <- dist/sum(dist)

  for(i in seq_len(maxtry)) {
    ## sample degree for each node among distribution
    degs <- sample(in1, size=n, replace=TRUE, prob=dist)
    ## if sum is not even, make it, by subtracting one from the last:
    if(sum(degs) %% 2 == 1)
      degs[which.max(degs)] <- degs[which.max(degs)] - 1L
    ## try: sometimes its not possible to construct graph for computed sequence
    g <- tryCatch(degree.sequence.game(degs, method="vl"), error = function(e) e)
    if(!inherits(g, "error"))
        break
  }
  if(i == maxtry && inherits(g, "error"))
      stop(gettextf("degree.sequence.game() did not succeed in maxtry=%d iterations",
                    maxtry), domain=NA)
  ## otherwise we are done
  simplify(g)
}


geoDAG <- function(n, r, dim, geo=TRUE, DAG, weighted,
                   wFUN = list(runif, min=0.1, max=1))
{
  stopifnot((n <- as.integer(n)) >= 2)
  points.dim <- matrix(runif(n*dim), dim, n)

  Q <- weights <- matrix(0, n, n)
  for(i in 1:(n-1)) {
    for(j in (i+1):n) {
      v <- points.dim[, i]-points.dim[, j]
      v <- v^2
      v <- sum(v)
      d <- sqrt(v)
      Q[i, j] <- d<=r

      c <- d/r
      c <- 1-c
      weights[i, j] <- (c*(1-0.1)+0.1)*Q[i, j]
    }
  }

  if(weighted) {
    if(geo) {
      Q <- weights
    } else {
      nrEdge <- sum(Q)
      Q[Q==1] <- do.call(wFUN[[1]], c(nrEdge, wFUN[-1]))
    }
  }

  if(!DAG) Q <- Q+t(Q)

  perm <- sample.int(n)
  Q <- Q[perm, perm]
  colnames(Q) <- rownames(Q) <- 1:n
  m2g(Q)
}


erDAG <- function(n, p)
{
  Q <- upper.tri(matrix(NA, n, n))
  Q[Q] <- rbinom((n-1)*n/2, 1, p)
  Q
}

generalHarmonic <- function(n, r) sum(1/(seq_len(n)^r))

findGamma <- function(n, d) {
  n1 <- n-1L
  uniroot(function(r) generalHarmonic(n1, r) / generalHarmonic(n1, r+1) - d,
          c(-10, 10))$root +1
}

##################################################
## unifDAG
##################################################
## A, B, a, sum, r, t: bigz

unifDAG <- function(n, weighted=FALSE, wFUN=list(runif, min=0.1, max=1)) {
  stopifnot(n>1)
  if (n > 100) stop("Use unifDAG only for n <= 100; for larger n use unifDAG.approx")

  ## step 1
  ## calculate numbers a_{n, k}, b_{n, k} and a_n up to N ##################
  ## is done offline #####################################


  ## step 2
  ## sample an integer between 1 and a_n ##########################
  r <- sampleZ2(.unifDagPreComp$a[n])

  ## step 3
  ## find vector K=c(k_1, ..., k_I) ##############################
  K <- findK.exact(n, r)

  ## step 4
  ## construct DAG (matrix Q) recursively ###########################
  Q <- sampleQ(n, K)

  if(weighted) {
    nrEdge <- sum(Q)
    if(!is.list(wFUN)) {wFUN <- list(wFUN)}
    Q[Q==1] <- do.call(wFUN[[1]], c(nrEdge, wFUN[-1]))
  }

  ## step 5
  ## permute matrix Q and convert to DAG #############################
  perm <- sample.int(n)
  as(Q[perm, perm], "graphNEL")

}


## find vector K=c(k_1, ..., k_I) ##############################
findK.exact <- function(n, r)
{
  K <- rep(0, n) # vector of k_1, ..., k_I
  k <- 1
  while(r>.unifDagPreComp$A[n, k]) {
    r <- r - .unifDagPreComp$A[n, k]
    k <- k+1
  }
  i <- 1
  K[i] <- k
  r <- as.bigz(as.bigq(r, chooseZ(n, k)))+1      #+1: should round to ceil
  m <- n-k
  while(m>0) {
    s <- 1
    t <- (2^k-1)^s * 2^as.bigz(k*(m-s)) * .unifDagPreComp$A[m, s]
    while(r>t) {
      r <- r-t
      s <- s+1
      if(m>=s) {t <- (2^k-1)^as.bigz(s) * 2^as.bigz(k*(m-s)) * .unifDagPreComp$A[m, s]}
      else {t <- r+1}
    }
    if(m>=s) {
      rn.z <- chooseZ(m, s) * (2^k-1)^as.bigz(s) * 2^as.bigz(k*(m-s))
      r.q <-  as.bigq(r, rn.z)
      r <- as.bigz(r.q)  + 1
      nn <- m
      k <- s
      m <- nn-k
      i <- i+1
      K[i] <- k}
    else {
      nn <- m
      k <- min(s, m)
      m <- nn-k
      i <- i+1
      K[i] <- k

    }
  }
  ## I <- i

  K[K!=0]
}

##' @title Sample Uniformly a Large (bigz) Integer
##' @param n a bigz (large) integer
##' @return a random large integer (class \code{"bigz"}) <= n
sampleZ2 <- function(n) {
### numbits <- as.integer(log2(n))+1
  numbits <- as.integer(log2(n-1))+1L
  repeat {
    r.bit <- rbinom(numbits, 1, prob=1/2) # from {0, 1}
    r <- as.bigz(paste0("0b", paste0(r.bit, collapse="")))
    if (r < n)
        return(r + 1)
  }
}


##################################################
## unifDAG.approx
##################################################
unifDAG.approx <- function(n, n.exact = 20, weighted=FALSE,
                           wFUN=list(runif, min=0.1, max=1)) {
  stopifnot(n>1)
  if (n < n.exact) stop("unifDAG.approx: n needs to be at least as big as n.exact!")

  ## step 1&2
  ## calculate numbers a_{n, k}, b_{n, k} and a_n up to N ##################
  ## calculate numbers A_k, B_{s|k} up to N.inf and accuracy #################
  ## is done offline #####################################

  ## step 3
  ## find approx-vector K=c(k_1, ..., k_I) #########################
  K <- findK.approx(n, n.exact)

  ## step 4
  ## construct DAG (matrix Q) recursively ###########################
  Q <- sampleQ(n, K)

  if(weighted) {
    nrEdge <- sum(Q)
    if(!is.list(wFUN)) {wFUN <- list(wFUN)}
    Q[Q==1] <- do.call(wFUN[[1]], c(nrEdge, wFUN[-1]))
  }

  ## step 5
  ## permute matrix Q and convert to DAG #############################
  perm <- sample.int(n)
  as(Q[perm, perm], "graphNEL")

}

## find vector K=c(k_1, ..., k_I) ##############################
findK.approx <- function(n, n.exact)
{
  M <- n
  K1 <- rep(0, n-n.exact)
  i <- 1
  K1[i]  <- sampleZ.cum.vec(.unifDagPreComp$Ak)
  M <- M-K1[i]
  i <- i+1
  while(M>n.exact) {
    K1[i] <- sampleZ.cum.vec(.unifDagPreComp$Bsk[, K1[i-1]])
    M <- M-K1[i]
    i <- i+1
  }
  if(M<n.exact) {
    M <- M+K1[i-1]
    K1[i-1] <- 0
  }
  K1 <- K1[K1!=0]

  K2 <- if(n.exact>=1) {
    ## direct enumeration method with n.exact
    r <- sampleZ2(.unifDagPreComp$a[M])
    findK.exact(M, r)
  }
  else
    0

  K <- c(K1, K2)
  K[K!=0]
}

sampleZ.cum.vec <- function(c) {
  ## c:: bigz-vector; c[i]=numbers of occurance of item i, returns random index, proportional to the numbers in c
  ind <- which(c!=0)
  s <- cumsum(c[ind])
  n <- length(ind)
  r <- sampleZ2(s[n])-1    ## since we want in [1:s[n]]
  ## linear search, since only small c is expected
  i <- 1
  while(s[i]<=r) {
    i <- i+1
  }
  ind[i]
}

##################################################
## Precompute data:
## List unifDagPreComp with elements
## Notation according to: Uniform random generation of large acyclic digraphs
## - A: a_{n, k}
## - B: b_{n, k}
## - a: a_n
## - Ak: A_k
## - Bsk: B_{s|k}
##################################################
if (FALSE) {
    library(gmp)
    setwd("/u/kalischm/research/packages/pcalg/pkg/R")
    source("genRandDAG.R")

    ## Exact
    load("/u/kalischm/research/packages/unifDAGs/tables100.RData")
    resExact <- generate.tables(100)
    ## identical(resExact[[1]], A) ## TRUE
    ## identical(resExact[[2]], B) ## TRUE
    ## identical(resExact[[3]], a) ## TRUE

    ## Approx
    load("/u/kalischm/research/packages/unifDAGs/tables_approx100_20.RData")
    resApprox <- approxK(N.inf=100, accuracy=20)
    ## identical(resApprox[[1]], Ak) ## TRUE
    ## identical(resApprox[[2]], Bsk) ## TRUE

    .unifDagPreComp <- list(A = resExact[[1]], B = resExact[[2]],
                           a = resExact[[3]],
                           Ak = resApprox[[1]], Bsk = resApprox[[2]])
    save(.unifDagPreComp,
         file = "/u/kalischm/research/packages/pcalg/pkg/sysdata.rda")
}

## calculate numbers a_{n, k}, b_{n, k} and a_n up to N ##################
## can be done offline ###################################
generate.tables <- function(N, dir=getwd(), verbose=TRUE) {

  A <- as.bigz(matrix(0, N, N))  # a_{n, k}
  B <- as.bigz(matrix(0, N, N))  # b_{n, k}
  a <- as.bigz(rep(0, N))       # a_n

  A[1, 1] <- B[1, 1] <- a[1] <- 1
  for(nn in 2:N) {
    if(verbose) cat("\n N: ", nn, " K: ")
    for(k in 1:(nn-1)) {
      if(verbose) cat(" ", k)
      sum <- as.bigz(0)
      for(s in 1:(nn-k)) {
        sum <- sum + (2^k-1)^as.bigz(s) * 2^as.bigz(k*(nn-k-s)) * A[nn-k, s]
      }
      B[nn, k] <-  sum
      A[nn, k] <- chooseZ(nn, k)*B[nn, k]
    }
    A[nn, nn] <- B[nn, nn] <- 1
    a[nn] <- sum(A[nn, 1:nn])
  }
  ## save(A, B, a, file=paste0(dir, "/tables", N, ".RData"))
  ## cat("\nTables saved in: ", paste0(dir, "/tables", N, ".RData"))
  list(A, B, a)
}

## construct A_k and B_{s|k} ###############################
approx.Ak <- function(N.inf=100, accuracy=20) {
  Ak <- as.bigz(rep(0, N.inf))  # A_k=lim_{n->oo}(A_{n, k}/a_n)

  acc <- 10^as.bigz(accuracy)
  for(k in 1:N.inf) {
    Ak[k] <- as.bigz(.unifDagPreComp$A[N.inf, k] * acc / .unifDagPreComp$a[N.inf])
  }
  Ak[Ak!=0]
}


approx.Bsk <- function(Ak) {
  n.k <- length(Ak)

  Bsk <- as.bigq(matrix(0, n.k, n.k))
  for(kk in 1:n.k) {
    for(ss in 1:n.k) {
      Bsk[ss, kk] <- as.bigq((1-1/(2^kk))^ss) * as.bigq(Ak[ss])
    }
  }
  as.bigz(Bsk)
}


## need table exact
approxK <- function(N.inf=100, accuracy=20, dir=getwd()) {
  Ak <- approx.Ak(N.inf, accuracy)
  Bsk <- approx.Bsk(Ak)

  ## save(Ak, Bsk, file=paste0(dir, "/tables_approx", N.inf, "_", accuracy, ".RData"))
  ## cat("\nApprox-Tables saved in: ", paste0(dir, "/tables_approx", N.inf, "_", accuracy, ".RData"))
  list(Ak, Bsk)
}


## augment a_{n, k}, b_{n, k} and a_n up form N0 to N ####################
## can be done offline ###################################
## augment.tables <- function(A0, B0, a0, N0, N, dir=getwd(), verbose=FALSE) {
##   A <- as.bigz(matrix(0, N, N))  # a_{n, k}
##   B <- as.bigz(matrix(0, N, N))  # b_{n, k}
##   a <- as.bigz(rep(0, N))       # a_n
##   A[1:N0, 1:N0] <- A0[1:N0, 1:N0]
##   B[1:N0, 1:N0] <- B0[1:N0, 1:N0]
##   a[1:N0] <- a0[1:N0]
##
##   for(nn in ((N0+1):N)) {
##     if(verbose) cat("\n N: ", nn, " K: ")
##     for(k in 1:(nn-1)) {
##       if(verbose) cat(" ", k)
##       sum <- as.bigz(0)
##       for(s in 1:(nn-k)) {
##         sum <- sum + (2^k-1)^as.bigz(s) * 2^as.bigz(k*(nn-k-s)) * A[nn-k, s]
##       }
##       B[nn, k] <-  sum
##       A[nn, k] <- chooseZ(nn, k)*B[nn, k]
##     }
##     A[nn, nn] <- B[nn, nn] <- 1
##     a[nn] <- sum(A[nn, 1:nn])
##   }
##
##   save(A, B, a, file=paste0(dir, "/tables", N, ".RData"))
##   cat("\nAugmented Tables saved in: ", paste0(dir, "/tables", N, ".RData"))
##   list(A, B, a)
## }
