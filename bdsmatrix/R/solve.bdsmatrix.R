# Cholesky decompostion and solution
solve.bdsmatrix<- function(a, b, full=TRUE, tolerance=1e-10, ...) {
    if (class(a) != 'bdsmatrix') 
	    stop("First argument must be a bdsmatrix")
    if (a@offdiag !=0) return(solve(as.matrix(a), b, tolerance=tolerance))
    nblock <- length(a@blocksize)
    adim <- dim(a)

    if (missing(b)) {
        # The inverse of the Cholesky is sparse, but if rmat is not 0
        #   the inverse of the martrix as a whole is not
        # For df computations in a Cox model, however, it turns out that
        #   I only need the diagonal of the matrix anyway.
        if (length(a@rmat)==0 || full==FALSE) {
            # The C-code will do the inverse for us
            temp <- .C("gchol_bdsinv", as.integer(nblock),
                                   as.integer(a@blocksize),
	                           as.integer(a@Dim),
                                   dmat= as.double(a@blocks),
                                   rmat= as.double(a@rmat),
                                   flag= as.double(tolerance),
	                           as.integer(0))

            if (length(a@rmat) >0) {
                new("bdsmatrix", blocksize=a@blocksize,
                    blocks = temp$dmat, offdiag=0, 
                    rmat = matrix(temp$rmat, nrow=nrow(a@rmat)),
                    Dim=a@Dim, Dimnames= a@Dimnames)
                }
            else {
                new("bdsmatrix", blocksize=a@blocksize,
                    blocks = temp$dmat, offdiag=0, 
                    Dim=a@Dim, Dimnames= a@Dimnames)
                }
            }
        else {
            # Get back the inverse of the cholesky from the C code
            #   and then multiply out the results ourselves (the C
            #   program doesn't have the memory space assigned to
            #   write out a full matrix).  The odds of a "not enough
            #   memory" message are high, if a is large.
            temp <- .C("gchol_bdsinv", as.integer(nblock),
                                   as.integer(a@blocksize),
	                           as.integer(a@Dim),
                                   dmat= as.double(a@blocks),
                                   rmat= as.double(a@rmat),
                                   flag= as.double(tolerance),
	                           as.integer(2))

            inv <- new('gchol.bdsmatrix', blocksize=a@blocksize, 
                       blocks=temp$dmat, 
                       rmat=matrix(temp$rmat, ncol=ncol(a@rmat)),
                       Dim=a@Dim, rank= as.integer(temp$flag),
                       Dimnames=a@Dimnames)
            dd <- diag(inv)
            rtemp <- as.matrix(inv)  #This may well complain about "too big"
            t(rtemp) %*% (dd* rtemp)
            }
	}
    
    else {
        #
        # If the rhs is a vector, save a little time by doing the decomp
        #  and the backsolve in a single .C call
        #
	if (length(b) == adim[1]) {
	    .C("gchol_bdssolve",as.integer(nblock),
	       		        as.integer(a@blocksize),
	       		        as.integer(a@Dim),
	       		        block = as.double(a@blocks),
	       		        rmat= as.double(a@rmat),
	       		        as.double(tolerance),
	       		        beta= as.double(b),
	                        flag=as.integer(0))$beta
	    }
	else {
            # The rhs is a matrix.  
            # In this case, it's faster to do the decomp once, and then
            #  solve against it multiple times.
            #
            if (!is.matrix(b) || nrow(b) != adim[1]) 
                stop("number or rows of b must equal number of columns of a")
            else solve(gchol(a, tolerance=tolerance), b)
            }
        }
    }
