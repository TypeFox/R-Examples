
##' Calculation of e^A with the Scaling & Squaring Method with Balancing
##' according to Higham (2008)
##'
##' R-Implementation of Higham's Algorithm from the Book (2008)
##' "Functions of Matrices - Theory and Computation", Chapter 10, Algorithm 10.20
##' Step 0:    Balancing
##' Step 1:    Scaling
##' Step 2:    Padé-Approximation
##' Step 3:    Squaring
##' Step 4:    Reverse Balancing
##'
##' @title Matrix Exponential with Scaling & Squaring and Balancing
##' @param A nxn Matrix
##' @param balancing logical indicating if balancing (step 0) should be applied
##' @return e^A Matrixeponential; nxn Matrix
##' @author Martin Maechler
expm.Higham08 <- function(A, balancing=TRUE)
{
    ## Check if A is square
    d <- dim(A)
    if(length(d) != 2 || d[1] != d[2]) stop("'A' must be a square matrix")
    n <- d[1]

    if (n <= 1) return(exp(A))

    ## else  n >= 2 ... non-trivial case : -------------

    ##---------STEP 0: BALANCING------------------------------------------------
    ## if balancing is asked for, balance the matrix A

    if (balancing) {
	baP <- balance(A,     "P")# -> error for non-classical matrix  -- "FIXME": balance()
	baS <- balance(baP$z, "S")
        A <- baS$z
    }


    ##--------STEP 1 and STEP 2 SCALING & PADÉ APPROXIMATION--------------------

    ## Informations about the given matrix
    nA <- Matrix::norm(A, "1")

    ## try to remain in the same matrix class system:
    I <- if(is(A,"Matrix")) Diagonal(n) else diag(n)

    ## If the norm is small enough, use the Padé-Approximation (PA) directly
    if (nA <= 2.1) {

    	t <- c(0.015, 0.25, 0.95, 2.1)
    	## the minimal m for the PA :
    	l <- which.max(nA <= t)

    	## Calculate PA
    	C <- rbind(c(120,60,12,1,0,0,0,0,0,0),
    		   c(30240,15120,3360,420,30,1,0,0,0,0),
    		   c(17297280,8648640,1995840,277200,25200,1512,56,1,0,0),
    		   c(17643225600,8821612800,2075673600,302702400,30270240,
    		     2162160,110880,3960,90,1))
        A2 <- A %*% A
    	P <- I
    	U <- C[l,2]*I
    	V <- C[l,1]*I

        for (k in 1:l) {
    	    P <- P %*% A2
    	    U <- U + C[l,(2*k)+2]*P
    	    V <- V + C[l,(2*k)+1]*P
    	}
    	U <- A %*% U
    	X <- solve(V-U,V+U)
    }

    ## Else, check if norm of A is small enough for m=13.
    ## If not, scale the matrix
    else {
        s <- log2(nA/5.4)
        B <- A
        ## Scaling
        if (s > 0) {
            s <- ceiling(s)
            B <- B/(2^s)
        }


	## Calculate PA

	c. <- c(64764752532480000,32382376266240000,7771770303897600,
                1187353796428800, 129060195264000,10559470521600, 670442572800,
                33522128640, 1323241920, 40840800,960960,16380, 182,1)

	B2 <- B %*% B
	B4 <- B2 %*% B2
	B6 <- B2 %*% B4

	U <- B %*% (B6 %*% (c.[14]*B6 + c.[12]*B4 + c.[10]*B2) +
		    c.[8]*B6 + c.[6]*B4 + c.[4]*B2 + c.[2]*I)
	V <- B6 %*% (c.[13]*B6 + c.[11]*B4 + c.[9]*B2) +
	    c.[7]*B6 + c.[5]*B4 + c.[3]*B2 + c.[1]*I

	X <- solve(V-U,V+U)

        ##---------------STEP 3 SQUARING----------------------------------------------
        if (s > 0) for (t in 1:s) X <- X %*% X
    }

    ##-----------------STEP 4 REVERSE BALANCING---------------------------------
    if (balancing) { ##	 reverse the balancing

	d <- baS$scale
	X <- X * (d * rep(1/d, each = n))

	## apply inverse permutation (of rows and columns):
	pp <- as.integer(baP$scale)
	if(baP$i1 > 1) {             ## The lower part
	    for(i in (baP$i1-1):1) { # 'p1' in *reverse* order
		tt <- X[,i]; X[,i] <- X[,pp[i]]; X[,pp[i]] <- tt
		tt <- X[i,]; X[i,] <- X[pp[i],]; X[pp[i],] <- tt
	    }
	}

	if(baP$i2 < n) {             ## The upper part
	    for(i in (baP$i2+1):n) { # 'p2' in *forward* order
		## swap	 i <-> pp[i]   both rows and columns
		tt <- X[,i]; X[,i] <- X[,pp[i]]; X[,pp[i]] <- tt
		tt <- X[i,]; X[i,] <- X[pp[i],]; X[pp[i],] <- tt
	    }
	}
    }

    X
}



##' Matrix Exponential -- using Al-Mohy and Higham (2009)'s algorithm
##'	 --> ../src/matexp_MH09.c
##' @param x square matrix
##' @param p the order of the Pade' approximation, 1 <= p <= 13.  The
##' default, 6, is what \file{expokit} uses.
expm.AlMoHi09 <- function(x, p = 6)
{
    d <- dim(x)
    if(length(d) != 2 || d[1] != d[2]) stop("'x' must be a square matrix")
    stopifnot(length(p <- as.integer(p)) == 1)
    if (p < 1 || p > 13)
        stop("Pade approximation order 'p' must be between 1 and 13.")

    .Call(R_matexp_MH09, x, p)
}
