setModel = function(x, model=NULL, chrom=NULL, penetrances=NULL, dfreq=NULL) {
	stopifnot(inherits(x,"linkdat"))
	if (!is.null(chrom)) stopifnot(is.character(chrom))
	if (!is.null(penetrances)) stopifnot(is.numeric(unlist(penetrances)), max(unlist(penetrances))<=1, min(unlist(penetrances))>=0)
	if (!is.null(dfreq)) stopifnot(is.numeric(dfreq), length(dfreq)==1, dfreq>=0, dfreq<=1)
	
	if (is.numeric(model)) {
		stopifnot(model %in% 1:4)
		model = switch(model, 
			list(chrom="AUTOSOMAL", penetrances=c(0,1,1), dfreq=1e-5), #aut dom
			list(chrom="AUTOSOMAL", penetrances=c(0,0,1), dfreq=1e-5), #aut rec
			list(chrom="X", penetrances=list(male=c(0,1), female=c(0,1,1)), dfreq=1e-5), #x-linked dom
			list(chrom="X", penetrances=list(male=c(0,1), female=c(0,0,1)), dfreq=1e-5) #x-linked rec
		)
	}
	if (is.null(model) && !is.null(x$model)) #if no model is given, but x already has one, use this as template
		model = x$model  
	hasmodel = !is.null(model)

	if (is.null(chrom)) chrom = ifelse(hasmodel, model$chrom, "AUTOSOMAL") else chrom <- match.arg(toupper(chrom), c("AUTOSOMAL","X"))
	if (is.null(dfreq)) dfreq = ifelse(hasmodel, model$dfreq, 1e-5)
	if (is.null(penetrances)) if (hasmodel) penetrances = model$penetrances else stop("No penetrance values given.")
	else {	switch(chrom,
		AUTOSOMAL = {
			if(is.character(penetrances)) {
				mod <- match.arg(tolower(penetrances), c("dominant", "recessive"))
				penetrances <- switch(mod, dominant=c(0,1,1), recessive=c(0,0,1))
			}
			if(length(penetrances)!=3) stop("For autosomal models, the penetrance parameter must be a vector of the form: c(f_0, f_1, f_2).")
			names(penetrances) <- c("f0","f1","f2")
		}, X = {
			if(is.character(penetrances)) {
				mod <- match.arg(tolower(penetrances), c("dominant", "recessive"))
				penetrances <- switch(mod, dominant=list(c(0,1), c(0,1,1)), recessive=list(c(0,1), c(0,0,1)))
			}
			if(any(length(penetrances)!=2, length(penetrances[[1]])!=2, length(penetrances[[2]])!=3)) stop("For X-linked models, the penetrance parameter must be a list of the form: list(c(f0_m, f1_m), c(f0_f, f1_f, f2_f)).")
			names(penetrances) <- c("male", "female");	names(penetrances$male) <- c("f0_m", "f1_m");	names(penetrances$female) <- c("f0_f", "f1_f", "f2_f")
		} )
	}

	#collecting the model information
	x$model = structure(list(chrom=chrom, penetrances=penetrances, dfreq=dfreq), class='linkdat.model')
	return(invisible(x))
}
	# #If and only if nallel==2 the following is carried out to give x an additional entry (x$initial_probs), containing initial likelihoods of each individual.
	# d=dfreq
	# switch(chrom,
	# AUTOSOMAL = {
		# p <- penetrances[c(3,2,1,3,2,2,1,3,2,1)]  # P(aff | geno). Note that P(non-aff | geno) = 1-p
		# Pen <- matrix( c( rep.int(1,10), 1-p, p), ncol = 3, dimnames = list(c('AADD','AADN','AANN','ABDD','ABDN','ABND','ABNN','BBDD','BBDN','BBNN'), 1:3))
		# DisFreq <- c(d^2, 2*d*(1-d), (1-d)^2, d^2, d*(1-d), d*(1-d), (1-d)^2, d^2, 2*d*(1-d), (1-d)^2)
		

		# penlist <- Pen[, x$pedigree[, 'AFF']+1]
		# penlist[, x$founders] <- penlist[, x$founders] * DisFreq 
		# #a=afreq[1]; b=afreq[2]; Afreq <- rep(c(a^2, 2*a*b, b^2), c(3, 4, 3))
		# #penlist[, x$founders] <- penlist[, x$founders] * Afreq
	# }, 
	# X = {
		# pM <- penetrances$male[c(2,1,2,1)] #P(aff | geno) for males
		# Pen_M <- matrix( c(rep.int(1,4), 1-pM, pM), ncol = 3,dimnames = list(c('AD','AN','BD','BN'), 1:3))

		# pF <- penetrances$female[c(3,2,1,3,2,2,1,3,2,1)] #P(aff | geno) for females
		# Pen_F <- matrix( c(rep.int(1,10), 1-pF, pF), ncol = 3,dimnames = list(c('AADD','AADN','AANN','ABDD','ABDN','ABND','ABNN','BBDD','BBDN','BBNN'), 1:3))
		
		# PenX <- list(male=Pen_M, female=Pen_F)
		# DisFreqX <- list(male=c(d, 1-d, d, 1-d), female=c(d^2, 2*d*(1-d), (1-d)^2, d^2, d*(1-d), d*(1-d), (1-d)^2, d^2, 2*d*(1-d), (1-d)^2))
		
		# ped=x$pedigree; sex=ped[, 'SEX']
		# penlist <- lapply(1:x$nInd, function(i) PenX[[ sex[i] ]][ , ped[i, 'AFF']+1])
		# for (i in x$founders) 	penlist[[i]] <- penlist[[i]] * DisFreqX[[ sex[i] ]] 
		# #a=afreq[1]; b=afreq[2]; AfreqX <- list(male=c(a, a, b, b), female=rep(c(a^2, 2*a*b, b^2), c(3, 4, 3)))
		# #for (i in x$founders) 	penlist[[i]] <- penlist[[i]] * AfreqX[[ sex[i] ]]
	# } )
	# x$initial_probs <- penlist
	# invisible(x)
# }
