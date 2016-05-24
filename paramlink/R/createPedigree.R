nuclearPed <- function(noffs, sex) {
	stopifnot(noffs>0)
	if(missing(sex)) sex = rep.int(1, noffs)
	if(length(sex) < noffs) sex = rep(sex, length.out=noffs)
	p = cbind(ID=1:(2+noffs), FID=c(0, 0, rep.int(1, noffs)), MID=c(0, 0, rep.int(2, noffs)),
			SEX=c(1, 2, sex), AFF=1)
	linkdat(p, verbose=FALSE)
}

cousinPed <- function(degree) {
	stopifnot(degree>=0)
	if(degree==0) return(nuclearPed(noffs=2, sex=1:2))
	p = cbind(ID=1:4, FID=c(0,0,1,1), MID=c(0,0,2,2), SEX=c(1,2,1,1), AFF=1)
	for (n in 1:degree)
		p = rbind(p, c(4*n+1,0,0,2,1), c(4*n+2,0,0,2,1), c(4*n+3,4*n-1,4*n+1,1,1), c(4*n+4,4*n,4*n+2,1,1))
	p[nrow(p), 'SEX'] = 2
	linkdat(p, verbose=FALSE)
}

halfCousinPed <- function(degree) {
	stopifnot(degree>=0)
	if(degree==0) p = cbind(ID=1:5, FID=c(0,0,0,1,1), MID=c(0,0,0,2,3), SEX=c(1,2,2,1,2), AFF=1)
	else {
		p = cbind(ID=1:3, FID=c(0,0,0), MID=c(0,0,0), SEX=c(1,2,2), AFF=1)
		for (n in seq_len(degree))
			p = rbind(p, c(4*n,4*n-4,4*n-3,1,1), c(4*n+1,0,0,2,1), c(4*n+2,4*n-2,4*n-1,1,1), c(4*n+3,0,0,2,1)) #add 1 generation: son in line 1, his wife, son in line 2, his wife
		dd = degree+1
		p = rbind(p, c(4*dd,4*dd-4,4*dd-3,1,1), c(4*dd+1,4*dd-2,4*dd-1,2,1)) #last generation - one boy, one girl.
		p[4,c(2,3)] = c(1,2); p[6,c(2,3)] = c(1,3) 
	}
	linkdat(p, verbose=FALSE)
}

mergePed <- function(x, y) {
    if(!is.null(x$markerdata) || !is.null(y$markerdata)) stop("Merging is only supported for pedigrees without marker data")
    ids = intersect(x$orig.ids, y$orig.ids)
    if(length(ids)==0) stop("Merging impossible: No common IDs")
    del = list(x = numeric(), y = numeric())
    for (i in ids) {
        if(.getSex(x,i) != .getSex(y,i)) stop(paste("Gender mismatch for individual", i))
        parx = parents(x, i)
        pary = parents(y, i)
        if(length(pary) == 0) del$y = c(del$y, i)
        else if(length(parx) == 0) del$x = c(del$x, i)
        else if(all(parx == pary)) del$y = c(del$y, i)
        else stop(paste("Parent mismatch for individual", i))
    }
    xx = as.matrix(x)[!x$orig.ids %in% del$x, ]
    yy = as.matrix(y)[!y$orig.ids %in% del$y, ]
    z = rbind(xx, yy)
    z[,'FAMID'] = x$famid # in case y$famid != x$famid
    
    # reorder to put parents above children (necessary when using IBDsim).
    N = nrow(z)
    i = 1
    while(i < N) {
        maxpar = max(match(z[i, c('FID','MID')], z[, 'ID'], nomatch=0))
        if(maxpar > i) {
            z = z[c(seq_len(i-1), (i+1):maxpar, i, seq_len(N-maxpar)+maxpar), ]
        }
        else i = i+1
    }
    
    restore_linkdat(z, attrs = attributes(xx))
}
