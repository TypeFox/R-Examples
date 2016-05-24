beta.temp <- function(x, y, index.family="sorensen"){

	# test for a valid index
	index.family <- match.arg(index.family, c('jaccard','sorensen'))

	# test for pre-existing betapart objects - validates data in x if not
	if (! inherits(x, "betapart")){	
		x <- betapart.core(x)
	}

	# test for pre-existing betapart objects - validates data in x if not
	if (! inherits(y, "betapart")){	
		y <- betapart.core(y)
	}

	# check x and y match
	if(! identical(dim(x$data), dim(y$data)))
		stop('The two data matrices do not have the same dimensions.')
	if(! identical(rownames(x$data), rownames(y$data)))
		stop('The two data matrices do not have the same site names.')
	if(! identical(colnames(x$data), colnames(y$data)))
		stop('The two data matrices do not have the same species names .')


	# get differences between groups
	ai<-apply(x$data &  y$data, 1, sum)
	bi<-apply(x$data & ! y$data, 1, sum)
	ci<-apply(! x$data & y$data, 1, sum)

	switch(index.family, 
		sorensen = {
			beta.sor<- (bi+ci) / (2*ai+bi+ci)
			beta.sim<- pmin(bi,ci)/(ai+pmin(bi,ci))
			beta.sne<- beta.sor - beta.sim
			result<-data.frame(beta.sim, beta.sne, beta.sor)},
		jaccard  = {
			beta.jac<- (bi+ci) / (ai+bi+ci)
			beta.jtu<- 2*pmin(bi,ci)/(ai+(2*pmin(bi,ci)))
			beta.jne<- beta.jac - beta.jtu
			result<-data.frame(beta.jtu, beta.jne, beta.jac)})

	return(result)
}
