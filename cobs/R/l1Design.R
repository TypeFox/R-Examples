####  Create B-Spline Design matrices for Sparse COBS :
####
####  l1.design2()  -- L_1        norm penalty ( <==> degree = 1 : order = 2)
####  loo.design2() -- L_Infinity norm penalty ( <==> degree = 2 : order = 3)
#### --------------------------------------------------------

### -->  For the moment, have  loo.design2()  moved to
### ./looDesign.R
###   ~~~~~~~~~~~   keep the two  "in parallel" !

l1.design2 <- function(x, w, constraint, ptConstr, knots,
                       pw, nrq, nl1, neqc, niqc, nvar, lambda)
{
    degree <- 1
    ks <- degree + 1 # k_{spline}, the spline "order"

    ##=########################################################################
    ##
    ## Generate the pseudo design matrix for L1 penalty
    ##
    ##=########################################################################

    ## create the pseudo design
    ##
    stopifnot(is.character(constraint),
              length(constraint) == length(niqc),
              nvar == as.integer(nvar), length(nvar) == 1, nvar >= 1)
    nvar <- as.integer(nvar)

    fieq <- FALSE ## 'fieq := we do Formulate InEQuality constraints'
    nk <- length(knots) + 2*(ks - 1) # = length(new.knots)
    ncoef <- nk - ks
    ## nobs == nrq + nl1 + neqc + niqc
    ox <- order(x)
    sortx <- x[ox]

    z1 <- .splBasis(ord = ks, knots, ncoef, xo = sortx)
    i1 <- rep.int(ox, rep.int(ks, nrq))
    j1 <- c(outer(1:ks, z1$offsets, "+"))
    Xeq <- as.matrix.csr(new("matrix.coo",
			     ra = c(t(t(z1$design)*rep(w[ox],ks))),
			     ia = i1, ja = j1,
			     dimension = as.integer(c(nrq, nvar))))
    nkm3 <- nk - 3
    nkm4 <- nk - 4
    niqc1 <- 0

    ##
    ## formulate the inequality constraints for the pseudo X  |--> Xieq
    ##
    if(lambda != 0 || !identical(constraint, "none")) {
	z2 <- .splBasis(ord = ks, knots, ncoef, xo = knots[1:nkm3],
			derivs = rep.int(1, nkm3))
	if(lambda != 0) {
	    ##
	    ## assign different weight to roughness penalty --- use 'pw' - FIXME!
            ## 'pw' is not used
	    ##
	    ra <- vector("numeric",(ks+1)*nl1)
	    if(nk > 4) { ## <==> nkm3 >= 2, nkm4 >= 1 :
		ra[seq(1,(ks+1)*nl1, by = 3)] <- -z2$design[1, 1:nkm4]
		ra[seq(3,(ks+1)*nl1, by = 3)] <-  z2$design[2, 2:nkm3]
		ra[seq(2,(ks+1)*nl1, by = 3)] <- (z2$design[1, 2:nkm3] -
						  z2$design[2, 1:nkm4])
	    }
	    ja <- c(outer(1:(ks+1),0:(nl1-1),"+"))
	    ia <- as.integer(seq(1, (ks+1)*(nl1+1), by = 3))
	    Xl1 <- new("matrix.csr", ra = ra * lambda, ja = ja, ia = ia,
		       dimension = as.integer(c(nl1, nvar)))
	    Xeq <- rbind(Xeq, Xl1)
	}

        ## TODO: start with empty Xieq and rbind(.) for every constraint
        ##  ---  BUT --- bug in 'SparseM' : cannot have a  '0 x 3' matrix.csr !
        if(FALSE) ## fails (BUG in 'SparseM') :
        Xieq <- new("matrix.csr", dimension = as.integer(c(0, nvar)),
                    ia = 1:1, ja=integer(0), ra = double(0))
        ## need to work with if(fieq) ... kludge:

	for(i.cnstr in seq(along=constraint)) {
	    constr <- constraint[i.cnstr]
	    niqc. <- niqc[i.cnstr]
	    if(constr == "increase" || constr == "decrease") {
		niqc4 <- nkm3
		ra <- c(if(constr == "increase") z2$design
				else -z2$design)
		ja <- as.integer(outer(1:ks, z2$offsets, "+"))
		ia <- as.integer(seq(1, ks*niqc4+1, by = ks))
		Xi <- new("matrix.csr", ra = ra, ja = ja, ia = ia,
			  dimension = as.integer(c(niqc4, nvar)))
		if(fieq)
		    Xieq <- rbind(Xieq, Xi)
		else {
		    Xieq <- Xi
		    fieq <- TRUE
		}
	    }
	    else if (constr == "periodic") {
		## *equality* constraint, coded as {+1, -1} INequalities

		## neqc3 <- 2
		z1.3 <- .splBasis(ord = ks, knots, ncoef,
				  xo = sortx[c(1, nrq)], derivs = c(1,1))
		ra <- c(z1.3$design)
		ja <- c(outer(1:ks, z1.3$offsets,"+"))
		ia <- as.integer(seq(1, length(z1.3$design)+1, by = ks))
		Xp <- new("matrix.csr", ra = ra, ja = ja, ia = ia,
			  dimension = as.integer(c(2, nvar)))
		Xp <- rbind(Xeq[ox[nrq],]- Xeq[ox[1],], Xp[2,]-Xp[1,])
		if(fieq)
		    Xieq <- rbind(Xieq, Xp, -Xp)
		else {
		    Xieq <- rbind(Xp,-Xp)
		    fieq <- TRUE
		}
	    }
	    else if(constr == "convex" || constr == "concave") {
		niqc4 <- nkm4
                if(nk > 4) { ## <==> niqc4 == nkm4 >= 1, nkm3 >= 2
                    sgn <- if(constr == "convex") +1 else -1
                    ra <- vector("numeric",(ks+1)*niqc4)
                    ra[seq(1,(ks+1)*niqc4, by = 3)] <- -z2$design[1,1:nkm4]*sgn
                    ra[seq(3,(ks+1)*niqc., by = 3)] <-  z2$design[2,2:nkm3]*sgn
                    ra[seq(2,(ks+1)*niqc., by = 3)] <- (z2$design[1,2:nkm3]-
                                           z2$design[2,1:nkm4])*sgn
                    ja <- c(outer(1:(ks+1),0:(niqc4-1),"+"))
                    ia <- as.integer(seq(1,(ks+1)*(niqc4+1), by = 3))
                    Xc <- new("matrix.csr", ra = ra, ja = ja, ia = ia,
                              dimension = as.integer(c(niqc., nvar)))
                    if(fieq)
                        Xieq <- rbind(Xieq, Xc)
                    else {
                        Xieq <- Xc
                        fieq <- TRUE
                    }
                }
		else
		    warning("too few knots ==> nk <= 4; could not add constraint '",constr,"'")
            }
        }## end{for}
    }

    if(ptConstr$n.smaller > 0) {
	o.smaller <- order(ptConstr$smaller[,2])
	smaller.o <- ptConstr$smaller[,2][o.smaller]
	z3.1 <- .splBasis(ord = ks, knots, ncoef, xo = smaller.o)
	ra <- -c(z3.1$design)
	ja <- c(outer(1:ks, z3.1$offsets,"+"))
	ia <- as.integer(seq(1, length(z3.1$design)+1, by = ks))
	Xp <- new("matrix.csr", ra = ra, ja = ja, ia = ia,
		  dimension = as.integer(c(ptConstr$n.smaller, nvar)))
	if(fieq)
	    Xieq <- rbind(Xieq, Xp)
	else {
	    Xieq <- Xp
	    fieq <- TRUE
	}
    }
    if(ptConstr$n.greater > 0) {
	o.greater <- order(ptConstr$greater[,2])
	greater.o <- ptConstr$greater[,2][o.greater]
	z3.2 <- .splBasis(ord = ks, knots, ncoef, xo = greater.o)
	ra <- c(z3.2$design)
	ja <- c(outer(1:ks, z3.2$offsets,"+"))
	ia <- as.integer(seq(1, length(z3.2$design)+1, by = ks))
	Xp <- new("matrix.csr", ra = ra, ja = ja, ia = ia,
		  dimension = as.integer(c(ptConstr$n.greater, nvar)))
	if(fieq)
	    Xieq <- rbind(Xieq, Xp)
	else {
	    Xieq <- Xp
	    fieq <- TRUE
	}
    }

    ##
    ## formulate the equality constraints
    ##
    ## Instead of addting to eq.constraints,
    ##      Xeq <- rbind(Xeq, Xp)
    ## we add it as + and -  to  Xieq :
    ##      Xieq. <- rbind(Xieq, Xp, -Xp)
    if(ptConstr$n.equal > 0) {
	o.equal <- order(ptConstr$equal[,2])
	equal.o <- ptConstr$equal[,2][o.equal]
	z1.1 <- .splBasis(ord = ks, knots, ncoef, xo = equal.o)
	ra <- c(z1.1$design)
	ja <- c(outer(1:ks, z1.1$offsets,"+"))
	ia <- as.integer(seq(1, length(z1.1$design)+1, by = ks))
	Xp <- new("matrix.csr", ra = ra, ja = ja, ia = ia,
		  dimension = as.integer(c(ptConstr$n.equal, nvar)))
	if(fieq)
	    Xieq <- rbind(Xieq, Xp,-Xp)
	else {
	    Xieq <- rbind(Xp,-Xp)
	    fieq <- TRUE
	}
    }

    if(ptConstr$n.gradient > 0) { ## gradient constraints for the pseudo X

	o.gradient <- order(ptConstr$gradient[,2])
	gradient.o <- ptConstr$gradient[,2][o.gradient]
	z1.2 <- .splBasis(ord = ks, knots, ncoef, xo = gradient.o,
			  derivs = rep.int(1, ptConstr$n.gradient))
	ra <- c(z1.2$design)
	ja <- c(outer(1:ks, z1.2$offsets,"+"))
	ia <- as.integer(seq(1, length(z1.2$design)+1, by = ks))
	Xp <- new("matrix.csr", ra = ra, ja = ja, ia = ia,
		  dimension = as.integer(c(ptConstr$n.gradient, nvar)))
	if(fieq)
	    Xieq <- rbind(Xieq, Xp,-Xp)
	else {
	    Xieq <- rbind(Xp,-Xp)
	    fieq <- TRUE
	}
    }
    if(fieq)
	return(list(Xeq = Xeq, Xieq = Xieq, fieq = fieq, niqc1 = niqc1))
    else
	return(list(Xeq = Xeq,              fieq = fieq, niqc1 = niqc1))
} ## l1.design
