### ===== File part of R package expm =====
###
### Function to compute the matrix exponential
###
###    exp(M) = sum(n = 0:Inf; M^n / n!),
###
### where M is an (n x n) matrix.
###

expm.s.Pade.s <- function(x, order, n=nrow(x)) {
    ## no checking here; this is not to be called by the user

    e <- ceiling(log2(max(rowSums(abs(x)))))
    s <- max(0, e+1)

    ## preconditions x :
    x <- x / (2^s)

    ## Pade approximation for exp(x)
    c <- .5
    D <- diag(1, n)
    E <- D + c*x
    D <- D - c*x
    X <- x
    p <- TRUE
    for(k in 2:order) {
	c <- c * (order-k+1) / (k*(2*order-k+1))
	X <- x %*% X # now  X = x ^ k
	cX <- c*X
	E <- E + cX
	D <- if(p) D + cX  else	 D - cX
	p <- !p
    }
    E <- solve(D, E)

    ## Undo the scaling by repeated squaring :
    for(k in seq_len(s))
	E <- E %*% E
    E
}


expm.methSparse <- c("Higham08", "R_Eigen", "R_Pade")
## "FIXME" -- keep this list up-to-date - test by setting  R_EXPM_NO_DENSE_COERCION

expm <- function(x, method = c("Higham08.b", "Higham08",
                    "AlMohy-Hi09",
		    "Ward77", "PadeRBS", "Pade", "Taylor", "PadeO", "TaylorO",
		    "R_Eigen", "R_Pade", "R_Ward77", "hybrid_Eigen_Ward"),
		 order = 8,
		 trySym = TRUE, tol = .Machine$double.eps,
		 do.sparseMsg = TRUE,
		 preconditioning = c("2bal", "1bal", "buggy"))
{
    ## some methods work for "matrix", "Matrix", or "mpfrMatrix", iff solve(.,.) worked:
    stopifnot(is.numeric(x) || (isM <- inherits(x, "dMatrix")) || inherits(x, "mpfrMatrix"))
    if(length(d <- dim(x)) != 2) stop("argument is not a matrix")
    if (d[1] != d[2]) stop("matrix not square")
    method <- match.arg(method)
    checkSparse <- !nzchar(Sys.getenv("R_EXPM_NO_DENSE_COERCION"))
    isM <- !is.numeric(x) && isM
    if(isM && checkSparse) { # i.e., a "dMatrix"
	if(!(method %in% expm.methSparse) && is(x, "sparseMatrix")) {
	    if(do.sparseMsg)
		message("coercing to dense matrix, as required by method ",
			dQuote(method))
	    x <- as(x, "denseMatrix")
	}
    }
    switch(method,
           "AlMohy-Hi09" = expm.AlMoHi09(x, p = order)
          ,
	   "Higham08.b" = expm.Higham08(x, balancing = TRUE)
	  ,
	   "Higham08"	= expm.Higham08(x, balancing = FALSE)
          ,
	   "Ward77" = {
	       ## AUTHORS: Christophe Dutang, Vincent Goulet at act ulaval ca
	       ##	 built on "Matrix" package, built on 'octave' code
	       ##	 Martin Maechler, for the preconditioning etc
               stopifnot(is.matrix(x))
	       switch(match.arg(preconditioning),
		      "2bal" = .Call(do_expm, x, "Ward77"),
		      "1bal" = .Call(do_expm, x, "Ward77_1"),
		      "buggy"= .Call(do_expm, x, "buggy_Ward77"),
		      stop("invalid 'preconditioning'"))
	   },
	   "R_Eigen" = {
	       ## matrix exponential using eigenvalues / spectral decomposition :
	       ## ==  Dubious Way 'Method 14' : is
	       ## good for 'symmetric' or 'orthogonal' (or other 'normal' : A'A = AA' ):

	       ## MM: improved from mexp2() with 'trySym' and isSymmetric()
	       isSym <- if(trySym) isSymmetric.matrix(x) else FALSE
	       z <- eigen(x, symmetric = isSym)
	       V <- z$vectors
	       Vi <- if(isSym) t(V) else solve(V)
	       Re(V %*% (    exp(z$values)   *	Vi)) ## ==
	       ##(V %*% diag(exp(z$values)) %*% Vi)
	   },
	   "hybrid_Eigen_Ward" = {
	       ## AUTHOR: Christophe Dutang
	       ## matrix exponential using eigenvalues / spectral decomposition and
	       ## Ward(1977) algorithm if x is numerically non diagonalisable
               stopifnot(is.matrix(x))
	       .Call(do_expm_eigen, x, tol)
	   },
	   "R_Pade"= { ## use scaling + Pade + squaring with R code:

	       ## matrix exponential using a scaling and squaring algorithm
	       ## with a Pade approximation
	       ## source code translated from expmdemo1.m in Matlab
	       ## by Stig Mortensen <sbm@imm.dtu.dk>,
	       ## prettified by MM -- works for "matrix" or "Matrix" matrices !

	       stopifnot(order >= 2)
	       expm.s.Pade.s(x, order, n=d[1])
	   },
	   "R_Ward77" = { ## R implementation of "Ward(1977)"
	       ## also works for "Matrix" matrices
	       stopifnot(order >= 2)
	       n <- d[1]
	       ## Preconditioning  Step 1: shift diagonal by average diagonal
	       trShift <- sum(d.x <- diag(x))
	       if(trShift) {
		   trShift <- trShift/n
		   diag(x) <- d.x - trShift
	       }

	       ## Preconditioning  Step 2: balancing with balance.
	       ##		   ------
	       ## For now, do like the octave implementation
	       ## TODO later:  use "B" (faster; better condition of result)
	       baP <- balance(x,     "P")
	       baS <- balance(baP$z, "S")

	       x <- expm.s.Pade.s(baS$z, order)
	       ##   -------------  scaling + Pade + squaring ------
	       ##   i.e., entails  Preconditioning Step 3 (and its reverse)

	       ## Reverse step 2: ------------------
	       ##
	       ## Step 2 b:  apply inverse scaling
	       d <- baS$scale
	       x <- x * (d * rep(1/d, each = n))
	       ##
	       ## Step 2 a:  apply inverse permutation (of rows and columns):
	       pp <- as.integer(baP$scale)
	       if(baP$i1 > 1) {		    ## The lower part
		   for(i in (baP$i1-1):1) { # 'p1' in *reverse* order
		       tt <- x[,i]; x[,i] <- x[,pp[i]]; x[,pp[i]] <- tt
		       tt <- x[i,]; x[i,] <- x[pp[i],]; x[pp[i],] <- tt
		   }
	       }
	       if(baP$i2 < n) {		    ## The upper part
		   for(i in (baP$i2+1):n) { # 'p2' in *forward* order
		       ## swap	 i <-> pp[i]   both rows and columns
		       tt <- x[,i]; x[,i] <- x[,pp[i]]; x[,pp[i]] <- tt
		       tt <- x[i,]; x[i,] <- x[pp[i],]; x[pp[i],] <- tt
		   }
	       }

	       ## reverse step 1 (diagonal shift)
	       if(trShift) {
		   exp(trShift) * x
	       }
	       else x
	   },
	   "PadeRBS" =  { ## the "expofit" method by  Roger B. Sidje (U.Queensland, AU)
	       stopifnot(is.matrix(x),
			 (order <- as.integer(order)) >= 1)
	       storage.mode(x) <- "double"
	       Fobj <- .Fortran(matexpRBS,
				order,		   # IDEG  1
				as.integer(d[1]),  # M	   2
				T = 1.,		   # T	   3
				H = x,		   # H	   4
				iflag = integer(1) # IFLAG 5
				)[c("H","iflag")]
	       if(Fobj[["iflag"]] < 0)
		   stop("Unable to determine matrix exponential")
	       Fobj[["H"]]
           }
	   , { ## the "mexp" methods by
	       ## AUTHORS: Marina Shapira and David Firth --------------
               stopifnot(is.matrix(x))
	       storage.mode(x) <- "double"
	       order <- as.integer(order)
	       ## MM:	a "silly"  way to code the method / order
	       ntaylor <- npade <- 0L
	       if (substr(method,1,4) == "Pade")
		   npade <- order else ntaylor <- order
	       res <- if(identical(grep("O$", method), 1L))
			  .Fortran(matrexpO,
			       X = x,
			       size = d[1],
			       ntaylor,
			       npade,
			       accuracy = double(1))[c("X", "accuracy")]
		      else
			  .Fortran(matrexp,
			       X = x,
			       size = d[1],
			       ntaylor,
			       npade,
			       accuracy = double(1))[c("X", "accuracy")]
	       structure(res$X, accuracy = res$accuracy)
	   })## end{switch}
}

