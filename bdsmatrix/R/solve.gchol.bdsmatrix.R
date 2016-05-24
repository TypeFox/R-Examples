# Backsolve or invert a gchol decompostion of a bds matrix
#  The "toler" arg to the C routines isn't used for this case, so
#  a dummy value of 0 has been inserted.  (Tolerance only is used in
#  the initial Cholesky decompostion).
# Assume that A is a bdsmatrix.  This routine mostly exists so that
#  solve(gchol(A), x) will give the same solution as solve(A,x).
#  Occasionally, the full=F argument may be needed as well.
#
solve.gchol.bdsmatrix<- function(a, b, full=TRUE, ...) {
    if (!inherits(a, 'gchol.bdsmatrix')) 
	    stop("First argument must be the gchol of a bdsmatrix")

    if (full) flag<-1  else flag <- 3
    nblock <- length(a@blocksize)
    if (length(a@rmat)==0) rmat <- 0.0  #dummy value to keep .C happy
    else rmat <- as.double(c(a@rmat))
    adim <- dim(a)

    if (missing(b)) {
	temp <- .C("gchol_bdsinv", as.integer(nblock),
                                   as.integer(a@blocksize),
		                   as.integer(a@Dim),
                                   dmat= as.double(a@blocks),
                                   rmat= rmat,
                                   as.double(0.0),
		                   as.integer(flag),
                                   copy=c(F,F,F,T,T,F,F))
	if (length(a@rmat) >0) {
            if (full)
                new('bdsmatrix',  blocksize=a@blocksize, blocks=temp$dmat,
                    rmat=matrix(temp$rmat, nrow=nrow(a@rmat)),
                    Dim=a@Dim, offdiag=0., Dimnames=a@Dimnames)
            else
                new('gchol.bdsmatrix', blocksize=a@blocksize, blocks=temp$dmat,
                     rmat=matrix(temp$rmat, nrow=nrow(a@rmat)),
                    Dim=a@Dim, rank=a@rank, Dimnames=a@Dimnames)
            }
	else {
            if (full)
                new('bdsmatrix', blocksize=a@blocksize, blocks=temp$dmat,
			Dim=a@Dim, offdiag=0., Dimnames=a@Dimnames)
            else
                new('gchol.bdsmatrix', blocksize=a@blocksize, blocks=temp$dmat,
			Dim=a@Dim, rank=a@rank, Dimnames=a@Dimnames)
            }
	}

    else {
	if (length(b) == adim[1]) {
	    .C("gchol_bdssolve",as.integer(nblock),
	       		        as.integer(a@blocksize),
	       		        as.integer(adim),
	       		        block = as.double(a@blocks),
	       		        rmat= rmat,
	       		        as.double(0.0),
	       		        beta= as.double(b),
	                        as.integer(flag),
                                copy=c(F,F,F,F,F,F,T,F))$beta
	    }
	else if (!is.matrix(b) || nrow(b) != adim[1]) 
	    stop("number or rows of b must equal number of columns of a")
	else {
	    temp <- b
	    for (i in 1:ncol(temp)) {
		temp[,i] <- .C("gchol_bdssolve",as.integer(nblock),
	       		        as.integer(a@blocksize),
	       		        as.integer(adim),
	       		        block = as.double(a@blocks),
	       		        rmat= rmat,
	       		        as.double(0.0),
	       		        beta= as.double(b[,i]),
	                        as.integer(flag),
                                copy=c(F,F,F,F,F,F,T,F))$beta
		}
	    temp
	    }
	}
    }
