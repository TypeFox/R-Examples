dgebalTst <- function(A) {

    ## Purpose: Consistency checking of	 dgebal()
    ## ----------------------------------------------------------------------
    ## Arguments: a square matrix
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, 20 Feb 2008 and on

    n <- dim(A)[1]
    ## do *the* three calls and look at result
    P <- dgebal(A, "P")

    doPerm <- function(A, pp, i1, i2) {
        stopifnot(length(pp) == n, dim(A) == c(n,n),
                  1 <= i1, i1 <= i2, i2 <= n)
        A. <- A
        if(i2 < n) { ## The upper part
            for(i in n:(i2+1)) {    # 'p2' in *reverse* order
                ## swap	 i <-> pp[i]   both rows and columns
                tt <- A.[,i]; A.[,i] <- A.[,pp[i]]; A.[,pp[i]] <- tt
                tt <- A.[i,]; A.[i,] <- A.[pp[i],]; A.[pp[i],] <- tt
            }
        }
        if(i1 > 1) { ## The lower part
            for(i in 1:(i1-1)) {    # 'p1' in *forward* order
                tt <- A.[,i]; A.[,i] <- A.[,pp[i]]; A.[,pp[i]] <- tt
                tt <- A.[i,]; A.[i,] <- A.[pp[i],]; A.[pp[i],] <- tt
            }
        }
        A.
    }

    checkPerm <- function(P, orig.A) {
        didPerm <- ((leftP  <- (i1 <- P$i1) != 1L) |
                    (rightP <- (i2 <- P$i2) != n))
        if(didPerm) { ## *had* permutation -- now check my idea about it
            pp <- as.integer(P$scale)
            ## Permute A to become P$z :
            A. <- doPerm(orig.A, pp = pp, i1=i1, i2=i2)
            stopifnot(isTRUE(all.equal(A., P$z, tolerance = 1e-15)))

            ## Now the reverse: Use pp[] and permute  A. "back to A":
            if(leftP) { ## The lower part
                for(i in (i1-1):1) {    # 'p1' in *reverse* order
                    tt <- A.[,i]; A.[,i] <- A.[,pp[i]]; A.[,pp[i]] <- tt
                    tt <- A.[i,]; A.[i,] <- A.[pp[i],]; A.[pp[i],] <- tt
                }
            }
            if(rightP) { ## The upper part
                for(i in (i2+1):n) {    # 'p2' in *forward* order
                    ## swap	 i <-> pp[i]   both rows and columns
                    tt <- A.[,i]; A.[,i] <- A.[,pp[i]]; A.[,pp[i]] <- tt
                    tt <- A.[i,]; A.[i,] <- A.[pp[i],]; A.[pp[i],] <- tt
                }
            }
            stopifnot(isTRUE(all.equal(A., orig.A, tolerance = 1e-15)))
        }
    }
    checkPerm(P, orig.A = A)

    S <- dgebal(P$z, "S")# "S" starting from result of "P"
    stopifnot(S$i1 == 1, S$i2 == n)

    ## Now check the scaling
    checkScal <- function (d, A1, A2) {
        stopifnot(length(d) == n, dim(A1) == dim(A2), dim(A2) == c(n,n))

        ## A.scaled <- diag(1/d, n) \%*\% A1 \%*\% diag(d, n)
        ## more efficiently:
        A.scaled <- A1 * (rep(d, each = n) / d)
        stopifnot(isTRUE(all.equal(A2, A.scaled, tolerance = 1e-15)))
        ## Check the reverse:
        S.rescaled <- A2 * (d * rep(1/d, each = n))
        stopifnot(isTRUE(all.equal(A1, S.rescaled, tolerance = 1e-15)))
    }
    checkScal(d = S$scale, A1 = P$z, A2 = S$z)

    B <- dgebal(A, "B")# "B" : B[oth]
    stopifnot(P$i1 == B$i1, P$i2 == B$i2)
    ## now check *both* permutation and scaling

    A.perm <- doPerm(A, pp = as.integer(B$scale), i1=B$i1, i2=B$i2)
    ## checkPerm(B, orig.A = A)

    dB <- B$scale
    dB[c(if(B$i1 > 1) 1:(B$i1-1),
         if(B$i2 < n) (B$i2+1):n)] <- 1
    checkScal(d = dB, A1 = A.perm, A2 = B$z)

    ## return
    list(P = P, S = S, B = B, Sz.eq.Bz = isTRUE(all.equal(S$z, B$z)))
}

m4. <- rbind(c(-1,-2, 0, 0),
             c( 0, 0,10,11),
             c( 0, 0,12, 0),
             c( 0,13, 0, 0))
str(b4. <- dgebalTst(m4.))

## better (?) example
(m <- matrix(c(0,-1,0,-2,10, rep(0,11)), 4,4))
str(ba <- dgebalTst(m))
## Hmm: here S$z  *differs*  from B$z
## ---  but at least, the scale[] and z[] returned seem ok


## a non-empty ``less-balanced'' example  ---

m4 <- matrix(outer(2^(0:7),c(-1,1)), 4,4)
m4[lower.tri(m4)] <- 0 #--> upper triangular ==> will have many permutations
## now permute it; so dgebal() will find the permutation
p <- c(4,2:1,3); m4 <- m4[p,p]
m4

str(dm4 <- dgebalTst(m4)) # much permutation!  i1 = i2 = 1 !
