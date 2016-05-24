##------OVERVIEW----------------------------------------------------------------

## Input:  A; nxn Matrix, no eigenvalues <=0, not singular
## Output: log(A); Matrixlogarithm; nxn Matrix


## Function for Calculation of log(A) with the Inverse Scaling&Squaring Method

## Step 0:    Schur Decompostion Tr
## Step 1:    Scaling (root of Tr)
## Step 2:    Padé-Approximation
## Step 3:    Squaring
## Step 4:    Reverse Schur Decomposition

## R-Implementation of Higham's Algorithm from the Book
## "Functions of Matrices - Theory and Computation", Chapter 11, Algorithm 11.9

##-------CODE-------------------------------------------------------------------

## The coefficients for the Padé-approximation can be computed at install time:

## r: exponents are in  (-51):(-56)
## p: exponents are in  c((-47):(-53), -56)


logm.H08.r <-
    rbind(c(5003999585967230*2^(-54), 8006399337547537*2^(-54), 5/18, 0,0,0,0),
          c(5640779706068081*2^(-51), 8899746432686114*2^(-53),
            8767290225458872*2^(-54), 6733946100265013*2^(-55), 0,0,0),
          c(5686538473148996*2^(-51), 4670441098084653*2^(-52),
            5124095576030447*2^(-53), 5604406634440294*2^(-54),
            8956332917077493*2^(-56), 0,0),
          c(5712804453675980*2^(-51), 4795663223967718*2^(-52),
            5535461316768070*2^(-53), 6805310445892841*2^(-54),
            7824302940658783*2^(-55), 6388318485698934*2^(-56), 0),
          c(5729264333934497*2^(-51), 4873628951352824*2^(-52),
            5788422587681293*2^(-53), 7529283295392226*2^(-54),
            4892742764696865*2^(-54), 5786545115272933*2^(-55),
            4786997716777457*2^(-56)))

logm.H08.p <-
    - rbind(c(7992072898328873*2^(-53), 1/2, 8121010851296995*2^(-56), 0,0,0,0),
            c(8107950463991866*2^(-49), 6823439817291852*2^(-51),
              6721885580294475*2^(-52), 4839623620596807*2^(-52), 0,0,0),
            c(6000309411699298*2^(-48),  4878981751356277*2^(-50), 2,
              5854649940415304*2^(-52), 4725262033344781*2^(-52),0,0),
            c(8336234321115872*2^(-48), 6646582649377394*2^(-50),
              5915042177386279*2^(-51), 7271968136730531*2^(-52),
              5422073417188307*2^(-52), 4660978705505908*2^(-52), 0),
            c(5530820008925390*2^(-47), 8712075454469181*2^(-50),
              7579841581383744*2^(-51), 4503599627370617*2^(-51),
              6406963985981958*2^(-52), 5171999978649488*2^(-52),
              4621190647118544*2^(-52)))


logm.Higham08 <- function(x) {
    ## work with "Matrix" too: x<-as.matrix(x)

    ##MM: No need to really check here; we get correct error msg later anyway
    ##    and don't need to compute det() here, in the good cases !
    ##    if (det(x) == 0) stop("'x' is singular")

    ##-------Step 0: Schur Decomposition-----------------------------------------

    ## Schur() checks for square matrix also:
    Sch.x <- Schur(Matrix(x, sparse=FALSE))
    ## FIXME 'sparse=FALSE' is workaround - good as long Matrix has no sparse Schur()
    ev <- Sch.x@EValues
    if(getOption("verbose") && any(abs(Arg(ev) - pi) < 1e-7))
## Let's see what works: temporarily *NOT* stop()ping :
	message(sprintf("'x' has negative real eigenvalues; maybe ok for %s", "logm()"))
    n <- Sch.x@Dim[1]
    Tr <- as.matrix(Sch.x@T)
    Q  <- as.matrix(Sch.x@Q)

    ##----- Step 1: [Inverse] Scaling -------------------------------------------
    I <- diag(n)
    thMax <- 0.264
    theta <- c(0.0162, 0.0539, 0.114, 0.187, thMax)
    p <- k <- 0 ; t.o <- -1
    ## NB: The following could loop forever, e.g., for  logm(Diagonal(x=1:0))
    repeat{
        t <- norm(Tr - I, "1") # norm(x, .) : currently x is coerced to dgeMatrix
	if(is.na(t)) {
	    warning(sprintf(ngettext(k,
				     "NA/NaN from %s after %d step.\n",
				     "NA/NaN from %s after %d steps.\n"),
			    " || Tr - I || ", k),
		    "The matrix logarithm may not exist for this matrix.")
	    return(array(t, dim=dim(Tr)))
	}
        if (t < thMax) {
            ## FIXME: use findInterval()
            j2 <- which.max( t <= theta)
            j1 <- which.max( (t/2) <= theta)
            if ((j2-j1 <= 1) || ((p <- p+1) == 2)) {
                m <- j2 ## m := order of the Padé-approximation
                break
            }
        } else if(k > 20 && abs(t.o - t) < 1e-7*t) {
            ##
	    warning(sprintf("Inverse scaling did not work (t = %g).\n", t),
		    "The matrix logarithm may not exist for this matrix.",
		    "Setting m = 3 arbitrarily.")
            m <- 3
            break
        }
        Tr <- rootS(Tr)##--> Matrix Square root of Jordan T
        ##    -----    [see below;  compare with ./sqrtm.R
        t.o <- t
        k <- k+1
    }
    if(getOption("verbose"))
	message(sprintf("logm.Higham08() -> (k, m) = (%d, %d)", k,m))

    ##------ Step 2: Padé-Approximation -----------------------------------------

    ## of order m :
    r.m <- logm.H08.r[m,]
    p.m <- logm.H08.p[m,]
    X <- 0
    Tr <- Tr-I
    for (s in 1:(m+2)) {
        X <- X + r.m[s]*solve(Tr - p.m[s]*I, Tr)
    }

    ##--- Step 3 & 4: Squaring & reverse Schur Decomposition -----------------

    2^k* Q %*% X %*% solve(Q)
}

### --- was  rootS.r -----------------------------------------------------------
###          ~~~~~~~

##------OVERVIEW----------------------------------------------------------------

## Input:  UT; nxn upper triangular block matrix (real Schur decomposition)
## Output: root of matrix UT, nxn upper triangular Matrix

## Function for calculation of UT^(1/2), which is used for the logarithm function

## Step 0:    Analyse block structure
## Step 1:    Calculate diagonal elements/blocks
## Step 2:    Calculate superdiagonal elements/blocks

## R-Implementation of Higham's Algorithm from the Book
## "Functions of Matrices - Theory and Computation", Chapter 6, Algorithm 6.7

rootS <- function(UT) {
    ## Generate Basic informations of Matrix UT
    stopifnot(length(d <- dim(UT)) == 2, is.numeric(d),
              (n <- d[1]) == d[2], n >= 1)
    ## FIXME : should work for "Matrix" too: not     S <- as.matrix(UT)
    S <- UT

    ##------- STEP 0: Analyse block structure ----------------------------------
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

    ##---------STEP 1: Calculate diagonal elements/blocks------------------------
    ## Calculate the root of the diagonal blocks of the Schur Decompostion S
    I <- diag(2)
    X <- matrix(0,n,n)
    for (j in seq_len(n-k)) {
        ij <- R.index[[j]]
	if (length(ij) == 1) {
	    ## Sij <- S[ij,ij]
	    ## if(Sij < 0)
	    ##	   ## FIXME(?) : in sqrtm(), we take *complex* sqrt() if needed :
	    ##	   ## -----  but afterwards  norm(Tr - I, "1") fails with complex
	    ##	   ## Sij <- complex(real = Sij, imaginary = 0)
	    ##	   stop("negative diagonal entry -- matrix square does not exist")
	    ## X[ij,ij] <- sqrt(Sij)
	    X[ij,ij] <- sqrt(S[ij,ij])
	}
        else {
            ## "FIXME"(better algorithm): only need largest eigen value
            ev1 <- eigen(S[ij,ij], only.values=TRUE)$values[1]
            r1 <- Re(sqrt(ev1)) ## sqrt(<complex>) ...
            X[ij,ij] <-
                r1*I + 1/(2*r1)*(S[ij,ij] - Re(ev1)*I)
        }
    }

### ___ FIXME __ code re-use: All the following is identical to 'STEP 3' in sqrtm()
###     -----    and almost all of STEP 1 above is == 'STEP 2' of sqrtm()

    ##---------STEP 2: Calculate superdiagonal elements/blocks-------------------

    ## Calculate the remaining, not-diagonal blocks
    if (n-k > 1) for (j in 2:(n-k)) {
        ij <- R.index[[j]]
        for (i in (j-1):1) {
            ii <- R.index[[i]]
            sumU <- 0

            ## Calculation for 1x1 Blocks
	    if (length(ij) == 1 & length(ii) == 1 ) {
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
			if (length(il) == 2) X[ii,il]%*%X[il,ij]
			else		     X[ii,il] * X[il,ij]
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
		X[ii,ij] <- solve(X[ii,ii]+X[ij,ij]*I,S[ii,ij]-sumU)
	    }
	    ## Calculation for	2x2 Blocks with special equation for solver
	    else if (length(ij) == 2 & length(ii) == 2 ) {
		if (j-i > 1) for (l in(i+1):(j-1)) {
		    il <- R.index[[l]]
		    sumU <- sumU + {
			if (length(il) == 2 ) X[ii,il] %*%   X[il,ij]
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
		X[ii,ij] <- solve(tUii+tUjj,as.vector(S[ii,ij]-sumU))
	    }
	}
    }
    X
}

