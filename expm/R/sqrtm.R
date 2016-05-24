####  Define sqrtm()    --- was Michael Stadelmann's  root.r
####         =======                                  ~~~~~~

##------OVERVIEW----------------------------------------------------------------

## Input:  A; nxn matrix, no eigenvalues <=0, not singular
## Output: root of matrix A, nxn Matrix


## Function for calculation of A^(1/2) with the real Schur decomposition

## Step 0:    real Schur decomposition T of A
## Step 1:    Aalyse block structure of T
## Step 2:    Calculate diagonal elements/blocks of T^(1/2)
## Step 3:    Calculate superdiagonal elements/blocks of T^(1/2)
## Step 4:    reverse Schur decompostion

## R-Implementation of Higham's Algorithm from the Book
## "Functions of Matrices - Theory and Computation", Chapter 6, Algorithm 6.7

sqrtm <- function(x) {
    ## Generate Basic informations of matrix x
    ## x <- as.matrix(x)
    d <- dim(x)
    if(length(d) != 2 || d[1] != d[2]) stop("'x' must be a quadratic matrix")

    ##MM: No need to really check here; we get correct error msg later anyway
    ##	  and don't need to compute det() here, in the good cases !
    ##	  if (det(x) == 0) stop("'x' is singular")
    n <- d[1]

    ##------- STEP 0: Schur Decomposition ---------------------------------------

    Sch.x <- Schur(Matrix(x)) ## <- {FIXME [Matrix]}
    ev <- Sch.x@EValues
    if(getOption("verbose") && any(abs(Arg(ev) - pi) < 1e-7))
        ## Let's see what works: temporarily *NOT* stop()ping :
        message(sprintf("'x' has negative real eigenvalues; maybe ok for %s", "sqrtm()"))

    S <- as.matrix(Sch.x@T)
    Q <- as.matrix(Sch.x@Q)

    ##---------STEP 1: Analyse block structure-----------------------------------

    ## Count 2x2 blocks (as Schur(x) is the real Schur Decompostion)
    J.has.2 <- S[cbind(2:n, 1:(n-1))] != 0
    k <- sum(J.has.2) ## := number of non-zero SUB-diagonals

    ## Generate Blockstructure and save it as R.index
    R.index <- vector("list",n-k)
    l <- 1
    i <- 1
    while(i < n) { ## i advances by 1 or 2, depending on 1- or 2- Jordan Block
	if (S[i+1,i] == 0) {
	    R.index[[l]] <- i
	}
	else {
	    R.index[[l]] <- (i:(i+1))
	    i <- i+1
	}
	i <- i+1
	l <- l+1
    }
    if (is.null(R.index[[n-k]])) { # needed; FIXME: should be able to "know"
        ##message(sprintf("R.index[n-k = %d]] is NULL, set to n=%d", n-k,n))
	R.index[[n-k]] <- n
    }

    ##---------STEP 2: Calculate diagonal elements/blocks------------------------
    ## Calculate the root of the diagonal blocks of the Schur Decompostion S
    I <- diag(2)
    X <- matrix(0,n,n)
    for (j in seq_len(n-k)) {
	ij <- R.index[[j]]
	if (length(ij) == 1) {
	    X[ij,ij] <- if((.s <- S[ij,ij]) < 0) sqrt(.s + 0i) else sqrt(.s)
	}
	else {
	    ev1 <- Sch.x@EValues[ij[1]]
	    r1 <- Re(sqrt(ev1)) ## sqrt(<complex>) ...
	    X[ij,ij] <- r1*I + 1/(2*r1)*(S[ij,ij] - Re(ev1)*I)
	}
    }
    ##---------STEP 3: Calculate superdiagonal elements/blocks-------------------

    ## Calculate the remaining, not-diagonal blocks
    if (n-k > 1) for (j in 2:(n-k)) {
	ij <- R.index[[j]]
	for (i in (j-1):1) {
	    ii <- R.index[[i]]
	    sumU <- 0

	    ## Calculation for 1x1 Blocks
	    if (length(ij) == 1 & length(ii) == 1) {
		if (j-i > 1) for (l in (i+1):(j-1)) {
		    il <- R.index[[l]]
		    sumU <- sumU + {
			if (length(il) == 2 ) X[ii,il]%*%X[il,ij]
			else		      X[ii,il] * X[il,ij]
		    }
		}
		X[ii,ij] <- solve(X[ii,ii]+X[ij,ij],S[ii,ij]-sumU)
	    }

	    ## Calculation for	1x2 Blocks
	    else if (length(ij) == 2 & length(ii) == 1 ) {
		if (j-i > 1) for (l in(i+1):(j-1)) {
		    il <- R.index[[l]]
		    sumU <- sumU + {
			if (length(il) == 2 ) X[ii,il]%*%X[il,ij]
			else		      X[ii,il] * X[il,ij]
		    }
		}
		X[ii,ij] <- solve(t(X[ii,ii]*I + X[ij,ij]),
                                  as.vector(S[ii,ij] - sumU))
	    }
	    ## Calculation for	2x1 Blocks
	    else if (length(ij) == 1 & length(ii) == 2 ) {
		if (j-i > 1) for (l in(i+1):(j-1)) {
		    il <- R.index[[l]]
		    sumU <- sumU + {
			if (length(il) == 2 ) X[ii,il]%*%X[il,ij]
			else		      X[ii,il] * X[il,ij]
		    }
		}
		X[ii,ij] <- solve(X[ii,ii]+X[ij,ij]*I, S[ii,ij]-sumU)
	    }
	    ## Calculation for	2x2 Blocks with special equation for solver
	    else if (length(ij) == 2 & length(ii) == 2 ) {
		if (j-i > 1) for (l in(i+1):(j-1)) {
		    il <- R.index[[l]]
		    sumU <- sumU + {
			if (length(il) == 2 ) X[ii,il] %*%  X[il,ij]
			else		      X[ii,il] %*% t(X[il,ij])

		    }
		}
		tUii <- matrix(0,4,4)
		tUii[1:2,1:2] <- X[ii,ii]
		tUii[3:4,3:4] <- X[ii,ii]
		tUjj <- matrix(0,4,4)
		tUjj[1:2,1:2] <- t(X[ij,ij])[1,1]*I
		tUjj[3:4,3:4] <- t(X[ij,ij])[2,2]*I
		tUjj[1:2,3:4] <- t(X[ij,ij])[1,2]*I
		tUjj[3:4,1:2] <- t(X[ij,ij])[2,1]*I
		X[ii,ij] <- solve(tUii+tUjj, as.vector(S[ii,ij]-sumU))
	    }
	}
    }

    ##------- STEP 4: Reverse the Schur Decomposition --------------------------
    ## Reverse the Schur Decomposition
    Q %*% X %*% solve(Q)
}

