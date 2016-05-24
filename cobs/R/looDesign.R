####  Create B-Spline Design matrices for Sparse COBS :
####
####  l1.design2()  -- L_1        norm penalty ( <==> degree = 1 : order = 2)
####  loo.design2() -- L_Infinity norm penalty ( <==> degree = 2 : order = 3)
#### --------------------------------------------------------

### -->  For the moment, have  l1.design2()  in
### ./l1Design.R
###   ~~~~~~~~~~    keep the two  "in parallel" !

loo.design2 <- function(x, w, constraint, ptConstr, knots,
                        pw, nrq, nl1, neqc, niqc, nvar, lambda)
{
    degree <- 2
    ks <- degree + 1 # k_{spline}, the spline "order"

    ##=########################################################################
    ##
    ## Generate the pseudo design matrix for L_oo penalty
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
			     ra = c(t(t(z1$design) * rep(w[ox], ks))),
			     ia = i1, ja = j1,
			     dimension = as.integer(c(nrq, nvar))))
    nd <- nk - 5
    ## nrql1 <- nrq + nl1
    ## nrleq <- nrql1 + neqc
    ## ==> nobs == nrleq + niqc

    ##
    ## formulate the inequality constraints for the pseudo X  |--> Xieq -- s''()
    ##
    if(lambda != 0) {
	niqc1 <- 2*nd
	z2 <- .splBasis(ord = ks, knots, ncoef, xo = knots[1:nd],
			derivs = rep.int(2, nd))
	Xl1 <- new("matrix.csr", ra = lambda,
		   ja = as.integer(nvar), ia = 1:2,
		   dimension = as.integer(c(1, nvar)))
	Xieq <- new("matrix.csr", ra = c(z2$design,-z2$design),
		    ja = c(rep.int(outer(1:ks, z2$offsets, "+"),2)),
		    ia = as.integer(seq(1,2*length(z2$design)+1, by = ks)),
		    dimension = as.integer(c(niqc1, nvar-1)))

        ## assign different weight to roughness penalty --- use 'pw' - FIXME!
	## 'pw' is not used
	Xeq <- rbind(Xeq, Xl1)
	Xieq <- cbind(Xieq, as.matrix.csr(rep.int(1, niqc1)))
	fieq <- TRUE
    } else niqc1 <- 0

    niqc4 <- niqc - (niqc1 + ptConstr$n.smaller + ptConstr$n.greater)

    for(i.cnstr in seq(along=constraint)) {
        constr <- constraint[i.cnstr]
        niqc4. <- niqc4[i.cnstr]

        if(constr == "convex" || constr == "concave") {
	    if(lambda == 0) ## z2 not yet above
		z2 <- .splBasis(ord = ks, knots, ncoef, xo = knots[1:niqc4.],
				derivs = rep.int(2, niqc4.))
	    ra <- c(if(constr == "convex")
			    z2$design else -z2$design)
	    ja <- c(outer(1:ks, z2$offsets,"+"))
	    ia <- as.integer(seq(1, length(z2$design)+1, by = ks))
	    Xc <- new("matrix.csr", ra = ra, ja = ja, ia = ia,
		      dimension = as.integer(c(niqc4., nvar)))
	    if(fieq)
		Xieq <- rbind(Xieq, Xc)
	    else {
		Xieq <- Xc
		fieq <- TRUE
	    }
	}
        else if(constr == "periodic") {
            ## *equality* constraint, coded as {+1, -1} INequalities

	    ## neqc3 <- 2
	    z1.3 <- .splBasis(ord = ks, knots, ncoef,
			    xo = sortx[c(1,nrq)], derivs = c(1,1))
	    ra <- c(z1.3$design)
	    ja <- c(outer(1:ks, z1.3$offsets,"+"))
	    ia <- as.integer(seq(1, length(z1.3$design)+1, by = ks))
	    Xp <- new("matrix.csr", ra = ra, ja = ja, ia = ia,
		      dimension = as.integer(c(2,nvar)))
	    Xp <- rbind(Xeq[ox[nrq],]- Xeq[ox[1],], Xp[2,]-Xp[1,])
	    if(fieq)
		Xieq <- rbind(Xieq, Xp,-Xp)
	    else {
		Xieq <- rbind(Xp,-Xp)
		fieq <- TRUE
	    }
	}
        else if(constr == "increase" || constr == "decrease") {
	    z3 <- .splBasis(ord = ks, knots, ncoef, xo = knots[1:niqc4.],
			    derivs = rep.int(1, niqc4.))
	    ra <- c(if(constr == "increase")
                    z3$design else -z3$design)
	    ja <- c(outer(1:ks, z3$offsets,"+"))
	    ia	<- as.integer(seq(1, length(z3$design)+1, by = ks))
	    X1 <- new("matrix.csr", ra = ra, ja = ja, ia = ia,
		      dimension = as.integer(c(niqc4., nvar)))
	    if(fieq)
		Xieq <- rbind(Xieq, X1)
	    else {
		Xieq <- X1
		fieq <- TRUE
	    }
	}

    } ## end{for}

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
} ## loo.design
