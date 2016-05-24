#### -------------------*- mode: R; kept-new-versions: 25; kept-old-versions: 20 -*-
#### Exponential Condition Number
#### ---------------------
#### Compute the Exponential Condition Number
#### ("1" and Frobenius-Norm) "exactly" and approximately.
####
#### All algorithms are based on the Fréchet derivative,
#### i.e., the expmCond() functions call expmFrechet()
#### for the calculation of the Fréchet derivative.

expmCond <- function(A, method = c("1.est", "F.est", "exact"),
                     expm = TRUE, abstol = 0.1, reltol = 1e-6,
                     give.exact = c("both", "1.norm", "F.norm"))
{
    ## Input:  A; nxn Matrix
    ## Output: list $ expmCondF: Exponentialconditionnumber Frobeniusnorm; scalar
    ##              $ expmCond1: Exponentialconditionnumber 1-Norm; scalar
    ##              $ expm:  e^A  Matrixexponential; nxn Matrix
    d <- dim(A)
    if(length(d) != 2 || d[1] != d[2] || d[1] <= 1)
        stop("'A' must be a square matrix of dimension at least 2")
    method <- match.arg(method)
    give.exact <- match.arg(give.exact)
    switch(method,
           "1.est" = .expmCond.1(A, expm = expm),
           "F.est" = .expmCond.F    (A, expm = expm, abstol=abstol, reltol=reltol),
           "exact" = .expmCond.X(A, expm = expm, give = give.exact),
           stop("invalid 'method'"))
}

### The former 4 files from Michi Stadelmann --- all in one file
##  byte  date        name
##  ---- ------------ ---------------
##  2006 Jan 30 12:12 expcond.r
##  2086 Jan 30 10:45 expcondest1.r
##  1782 Jan 30 10:45 expcondestfrob.r
##  4544 Jan 30 12:22 expm2frech.r

###------------------ expcond.r -------------------------------------------

## Function for *eXact* (slow!) calculation of the Exponentialconditionnumber
## ("1" and Frobenius-Norm).
## R-Implementation of Higham's Algorithm from the book
## "Functions of Matrices - Theory and Computation", chapter 3.4, algorithm 3.17

## Step 1:    Calculate Kroneckermatrix of L(A)
## Step 2:    Calculate Expentialconditionnumber ("1" & Frobenius-Norm)

.expmCond.X <- function(A, give= c("both", "1.norm", "F.norm"), expm = TRUE)
{
    ## Input:  A; nxn Matrix
    ## Output: list $ expmCondF: Exponentialconditionnumber Frobeniusnorm; scalar
    ##              $ expmCond1: Exponentialconditionnumber 1-Norm; scalar
    ##              $ expm:  e^A  Matrixexponential; nxn Matrix
    d <- dim(A)
    if(length(d) != 2 || d[1] != d[2] || d[1] <= 1)
	stop("'A' must be a square matrix of dimension at least 2")
    n <- d[1]
    ##---------STEP 1: Calculate Kroneckermatrix of L(A)------------------------
    K  <- matrix(0, n^2, n^2)
    E0 <- matrix(0, n,n)
    E.unit <- function(i,j) {
        ## Compute E_ij in R^{n x n} , the ij-th unit Matrix
        E <- E0
        E[i,j] <- 1
        E
    }

    give <- match.arg(give)
    jj <- 0
    for (j in  1:n) {
	for (i in 1:n) {
	    calc <- expmFrechet(A, E.unit(i,j), expm=(j == n) && (i == n))
	    K[, (jj <- jj + 1)] <- calc$Lexpm
	}
    }
    ##-------STEP 2 CALCULATE EXPONENTIALCONDITIONNUMBER ------------------------

    ## Frobenius-Norm
    do.F <- (give %in% c("F.norm", "both"))
    do.1 <- (give %in% c("1.norm", "both"))
    if(do.F)
        normk <- sqrt(max(eigen(crossprod(K))$values)) # crossprod(K) := K' K
    list(expmCondF = ## Frobenius Norm
         if(do.F) normk      * norm(A,"F") /  norm(calc$expm,"F"),
         expmCond1 = ## 1-Norm
         if(do.1) norm(K,"1")* norm(A,"1") / (norm(calc$expm,"1")*n),
         expm = if(expm) calc$expm)
}

###------------------ expcondest1.r ---------------------------------------

## Function for Estimation of the "1"-norm exponentialcondtionnumber based on
## the LAPACK marix norm estimator.

## R-Implementation of Higham's Algorithm from the book
## "Functions of Matrices - Theory and Computation", chapter 3.4, algorithm 3.21

## Step 1:    Estimate "1"-Norm of Kroneckermatrix K(A)
##            This step is based on the equation: K(A)vec(E)=vec(L(A,E))
## Step 2:    Calculate Expentialconditionnumber ("1"-Norm)

.expmCond.1 <- function(A, expm = TRUE)
{
    ## Input:  A; nxn Matrix
    ## Output: list $ expmCond: Exponentialconditionnumber "1"-Norm; scalar
    ##              $ expm:     e^A   Matrixexponential; nxn Matrix

    ##-------STEP 1 ESTIMATE "1"-NORM FROM THE KRONECKERMATRIX--------------

    ## Check if A is square
    d <- dim(A)
    if(length(d) != 2 || d[1] != d[2] || d[1] <= 1)
	stop("'A' must be a square matrix of dimension at least 2")
    n <- d[1]

    tA <- t(A)
    E    <- matrix(1/n^2, n,n)
    calc <- expmFrechet(A,E)
    V    <- calc$Lexpm
    G    <- sum(abs(V))
    Z    <- sign(V)
    X    <- expmFrechet(tA,Z, expm=FALSE)$Lexpm

    k <- 2
    E0 <- matrix(0, n,n)
    repeat { ## at most steps k = 2, 3, 4, 5
        j <- which.max(as.vector(abs(X)))
        Ej <- E0; Ej[j] <- 1
        V  <- expmFrechet(A,Ej, expm=FALSE)$Lexpm
        G  <- sum(abs(V))
        sV <- sign(V)
        if (identical(sV, Z) ||
            identical(sV,-Z)) break
        Z     <- sV
        X     <- expmFrechet(tA,Z, expm=FALSE)$Lexpm
        k     <- k+1
        if (k > 5 || max(abs(X)) == X[j])
            break
    }
    ## 'G' = gamma now is our desired lower bound

    ## Now, try another "lucky guess" and

    ## *increase* G if the guess *was* lucky :
    for (l in 1:(n^2)) { ## FIXME: vectorize this!
        X[l] <- (-1)^(l+1) * (1+(l-1)/(n^2-1))
    }
    X <- expmFrechet(A,X, expm=FALSE)$Lexpm
    G. <- 2*sum(abs(X))/(3*n^2)
    if (G < G.) {
        message("'lucky guess' was better")
        G <- G.
    }

    ##-------STEP 2 CALCULATE EXPONENTIALCONDITIONNUMBER------------------
    C1 <- G * norm(A,"1") / (norm(calc$expm,"1")*n)
    if(expm) list(condExpm = C1, expm = calc$expm) else C1
}

###------------------ expcondestfrob.r ------------------------------------

## Function for estimation of the frobenius-Norm exponentialcondtionnumber based
## on the powermethod-matrixnorm estimation.

## R-Implementation of Higham's Algorithm from the book
## "Functions of Matrices - Theory and Computation", chapter 3.4, algorithm 3.19

## Step 1:    Estimate "2"-Norm of Kroneckermatrix K(A)
##            This step is based on the equation: K(A)vec(E)=vec(L(A,E))
## Step 2:    Calculate Expentialconditionnumber (Frobenius-Norm)

.expmCond.F <- function(A, abstol = 0.1, reltol = 1e-6,
                        maxiter = 100, expm = TRUE)
{
    ## Input:  A; nxn Matrix
    ## Output: list C: C$expmCond: Exponentialconditionnumber Frobeniusnorm; scalar
    ##                C$expm:    e^A Matrixexponential; nxn Matrix

    ## Check if A is square
    d <- dim(A)
    if(length(d) != 2 || d[1] != d[2] || d[1] <= 1)
	stop("'A' must be a square matrix of dimension at least 2")
    n <- d[1]

    ##-------STEP 1 ESTIMATE 2-NORM OF KRONECKERMATRIX-------------------------------
    Z1 <- if(is(A,"Matrix")) Matrix(rnorm(n*n),n,n) else matrix(rnorm(n*n),n,n)
    tA <- t(A)
    calc <- expmFrechet(A,Z1)
    W1   <- calc$Lexpm
    Z1   <- expmFrechet(tA,W1, expm=FALSE)$Lexpm
    G2   <- norm(Z1,"F")/norm(W1,"F")

    it <- 0
    repeat {
        G1 <- G2
        W2 <- expmFrechet(A, Z1, expm=FALSE)$Lexpm
        Z2 <- expmFrechet(tA,W2, expm=FALSE)$Lexpm
        G2 <- norm(Z2,"F")/norm(W2,"F")
        Z1 <- Z2
        dG <- abs(G1-G2)
        it <- it+1
        if (it > maxiter || dG < abstol && dG < reltol*G2)
            break
    }
    if(it > maxiter)
	warning(sprintf("reached  maxiter = %d iterations; tolerances too small?",
			maxiter))

    ##-------STEP 2 CALCULATE EXPONENTIALCONDITIONNUMBER--------------------
    cF <- G2*norm(A,"F") / norm(calc$expm,"F")
    attr(cF, "iter") <- it
    if(expm) list(condExpm = cF, expm = calc$expm) else cF
}

###------------------ expm2frech.r ----------------------------------------------

## Calculation of e^A and the Exponential Frechet-Derivation L(A,E)
## with the Scaling & Squaring Method

## R-Implementation of Higham's Algorithm from the Article
## "Computing Fréchet Derivative of the Matrix Exponential, with an application
## to Condition Number Estimation", MIMS EPrint 2008.26, Algorithm 6.4

## Step 1:    Scaling (of A and E)
## Step 2:    Padé-Approximation of e^A and L(A,E)
## Step 3:    Squaring

expmFrechet <- function(A,E, method = c("SPS","blockEnlarge"), expm = TRUE)
{
    ## Input:  A; nxn Matrix
    ##         E; nxn Matrix
    ## Output: list X: X$expm; e^A Matrixeponential; nxn Matrix
    ##                X$Lexpm; Exponential-Frechet-Derivative L(A,E); nxn Matrix

    ## Check if A is square
    d <- dim(A)
    if(length(d) != 2 || d[1] != d[2]) stop("'A' must be a square matrix")
    stopifnot(is.matrix(E))
    if(!identical(d,dim(E))) stop("A and E need to have the same dimension")
    n <- d[1]

    if (n <= 1) {
        X <- exp(A)
        X2<- E*X
        return(if(expm) list(expm= X, Lexpm = X2) else list(Lexpm = X2))
    }
    ## else  n >= 2 ... non-trivial case : -------------

    method <- match.arg(method)

    switch(method,
           "SPS" = .expmFrechet2008.26(A,E, expm = expm)
           ,
           "blockEnlarge" = {
               ## From: Daniel Kressner  @ math ETH Zurich
               ## To:   Stadelmann Michael, Cc: Martin Maechler
               ## Subject: Frechet-Ableitung von f testen
               ## Date: Mon, 26 Jan 2009

               ## mir ist noch ein weiterer Weg zum Test Deines
               ## Algorithmus fuer die Frechet-Ableitung eingefallen.

               ## Berechnet man  f ([A E, 0 A])

               ## dann enthaelt der (1,2)-Block die Ableitung von f an
               ## der Stelle A in Richtung E (siehe Higham).
               OO <- array(0, dim=d)
               B <- rbind(cbind(A,  E),
                          cbind(OO, A)) ## stopifnot(dim(B) == 2*d)
               fB <- expm.Higham08(B)[1:n, ]
               L <- fB[ , n+ 1:n]
               if(expm) list(expm = fB[ , 1:n], Lexpm = L) else list(Lexpm = L)
           })
} ## expmFrechet


.expmFrechet2008.26 <- function(A, E, expm = TRUE)
{
    ## No error checking!   --> not to be called by the user!

    ## R-Implementation of Higham's Algorithm from the Article
    ## "Computing Fréchet Derivative of the Matrix Exponential, with an application
    ## to Condition Number Estimation", MIMS EPrint 2008.26, Algorithm 6.4

    ## Step 1:    Scaling (of A and E)
    ## Step 2:    Padé-Approximation of e^A and L(A,E)
    ## Step 3:    Squaring

    ##-----------STEP 1 & STEP 2: SCALING & PADÉ APPROXIMATION-------------------

    ## Informations about the given matrix
    nA  <- norm(A ,"1") ## == Matrix::norm

    n <- nrow(A)# == ncol(A) .. tested "in the caller"
    ## try to remain in the same matrix class system:
    I <- if(is(A,"Matrix")) Diagonal(n) else diag(n)

    ## If the norm is small enough, use directly the Padé-Approximation (PA)
    if (nA <= 1.78) {

        t <- c(0.0108,0.2,0.783,1.78)
       	## the minimal m for the PA :
        l <- which.max(nA <= t)

        ## Calculate PA for e^A and L(A,E)
        C <- rbind(c(120,60,12,1,0,0,0,0,0,0),
                   c(30240,15120,3360,420,30,1,0,0,0,0),
                   c(17297280,8648640,1995840,277200,25200,1512,56,1,0,0),
                   c(17643225600,8821612800,2075673600,302702400,30270240,
                     2162160,110880,3960,90,1)) [l , ] # only need l-th row

        P  <- I
        U  <- C[2]*I
        V  <- C[1]*I
        A2 <- A %*% A
        M2 <- A %*% E + E %*% A
        M  <- M2
        LU <- C[4]*M
        LV <- C[3]*M

        oC <- 2
        for (k in seq_len(l-1)) { ## oC == 2k
            ## PA e^A
            P  <- P %*% A2
            U  <- U+C[oC+ 2]*P
            V  <- V+C[oC+ 1]*P

            ## PA L(A,E)
            M  <- A2 %*% M + M2 %*% P
            LU <- LU + C[oC+ 4]*M
            LV <- LV + C[oC+ 3]*M
            oC <- oC + 2
        }
        ## PA e^A & L(A,E)
        P  <- P %*% A2
        U  <- U + C[oC+ 2]*P
        LU <- A %*% LU + E %*% U
        U  <- A %*% U
        V  <- V + C[oC+ 1]*P

        X  <- solve(V-U, V+U)
        X2 <- solve(V-U, LU+LV + (LU-LV)%*%X)
    }

    ## Else, check if norm of A is small enough for PA with m=13.
    ## If not, scale the matrix
    else {
        s  <- log2(nA/4.74)
        B  <- A
        D  <- E
        ## Scaling
        if (s > 0){
            s  <- ceiling(s)
            B  <- A/(2^s)
            D  <- D/(2^s)
        }
        C. <- c(64764752532480000,32382376266240000,7771770303897600,1187353796428800,
                129060195264000,10559470521600,670442572800,33522128640,1323241920,
                40840800,960960,16380,182,1)

        ## Calculate PA
        ## PA e^A
        B2  <- B%*%B
        B4  <- B2%*%B2
        B6  <- B2%*%B4
        W1  <- C.[14]*B6+ C.[12]*B4+ C.[10]*B2
        W2  <- C.[ 8]*B6+ C.[ 6]*B4+ C.[ 4]*B2+C.[2]*I
        Z1  <- C.[13]*B6+ C.[11]*B4+ C.[ 9]*B2
        Z2  <- C.[ 7]*B6+ C.[ 5]*B4+ C.[ 3]*B2+C.[1]*I
        W   <- B6%*%W1+W2
        U   <- B%*%W
        V   <- B6%*%Z1+Z2

        ## PA L(A,E)
        M2  <- B%*%D + D%*%B
        M4  <- B2%*%M2 + M2%*%B2
        M6  <- B4%*%M2 + M4%*%B2
        LW1 <- C.[14]*M6+ C.[12]*M4+ C.[10]*M2
        LW2 <- C.[ 8]*M6+ C.[ 6]*M4+ C.[ 4]*M2
        LZ1 <- C.[13]*M6+ C.[11]*M4+ C.[ 9]*M2
        LZ2 <- C.[ 7]*M6+ C.[ 5]*M4+ C.[ 3]*M2
        LW  <- B6%*%LW1 + M6%*%W1 + LW2
        LU  <- B%*%LW + D%*%W
        LV  <- B6%*%LZ1 + M6%*%Z1 + LZ2

        X   <- solve(V-U, V+U)
        X2  <- solve(V-U, LU+LV + (LU-LV)%*%X)

        ##----------STEP 3 SQUARING----------------------------------------------
        ## Squaring
        if (s > 0) for (t in seq_len(s)) {
            X2 <- X2 %*% X  +  X %*% X2
            if(expm || t != s)
                X  <- X %*% X
        }
    }
    if(expm) list(expm = X, Lexpm = X2) else list(Lexpm = X2)
} ## .expmFrechet2008.26
