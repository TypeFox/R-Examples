f.sim <- function(.prob, size, nall, .nloci, xchrom, .grid){
##
## SAMPLE FROM MULTINOMIAL DISTRIBUTION
##	.sim.rownum <- rMultinom(matrix(.prob, nrow = 1), m = n.cases)
	.sim.rownum <- suppressWarnings(sample(length(.prob), size = size, replace = T, prob = .prob))
	#
	## FIND ALLELES CORRESPONDING TO GRID ROW NUMBERS
	if(xchrom){
		.nrows <- dim(.grid)[1]/2
		.girls <- (.sim.rownum > .nrows)
		.tmp.boys <- .sim.rownum[!.girls]
		.tmp.girls <- .sim.rownum[.girls]
		.tmp.girls <- .tmp.girls - .nrows
		.alleles.boys <- f.pos.to.haplocomb(A = nall, pos = .tmp.boys, fam = "mfx")
		.alleles.girls <- f.pos.to.haplocomb(A = nall, pos = .tmp.girls, fam = "mfx")
		###.sex1 <- .sex[.sim.rownum]
		.sex1 <- c(rep(1, length(.tmp.boys)), rep(2, length(.tmp.girls)))
		###.alleles <- rbind(.alleles.boys, .alleles.girls)
	} else{
		.alleles <- f.pos.to.haplocomb(A = nall, pos = .sim.rownum)
	}
	#
	## ADD COLUMNS WITH CHILD GENOTYPES
	if(xchrom){
		.names <- dimnames(.alleles.boys)[[2]]
		.names <- matrix(.names, nrow = .nloci)
		.ind.boys <- c(1,2,3,3,2,2)
		.ind.girls <- c(1,2,3,3,2,3)
		.names.boys <- .names[, .ind.boys]
		.names.girls <- .names[, .ind.girls]
		.names.boys <- as.vector(t(.names.boys))
		.names.girls <- as.vector(t(.names.girls))
		#
		.alleles.boys <- .alleles.boys[,.names.boys]
		.alleles.girls <- .alleles.girls[,.names.girls]
		.alleles <- rbind(.alleles.boys, .alleles.girls)
	}else{
		.names <- dimnames(.alleles)[[2]]
		.names <- matrix(.names, nrow = .nloci)
		.ind <- c(1,2,3,4,2,4)
		.names <- .names[, .ind]
		.names <- as.vector(t(.names))
		#
		.alleles <- .alleles[,.names]
	}
	#
	## RANDOMIZE SEQUENCE OF ALLELES AND PREPARE FOR WRITING TO DISK
	.all.paste <- vector(dim(.alleles)[2]/2, mode = "list")
	for(j in seq(along = .all.paste)){
		.all.paste[[j]] <- f.rand.geno(.alleles[,2*j-1], .alleles[,2*j])
	}

	.all.paste <- as.data.frame(.all.paste)
	names(.all.paste) <- seq(along = .all.paste)

	if(xchrom) .ut <- cbind(sex = .sex1, .all.paste)
	else .ut <- .all.paste
	#
	return(.ut)	
}	
