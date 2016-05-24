.peelSNP <- function(dat, sub, SEX, chrom, TR.MATR) {
	probs = dat[["probs"]];   haps = dat[["haps"]];   fa_haps = haps[[sub[["father"]]]];   mo_haps = haps[[sub[["mother"]]]];   fa_len = length(fa_haps);   mo_len = length(mo_haps)
	likel = probs[[sub[["father"]]]] %*% t.default(probs[[sub[["mother"]]]])
	piv = sub[["pivot"]]; pivtype = sub[["pivtype"]]
	offs = sub[['offspring']]; nonpiv.offs = offs[offs != piv]
	switch(chrom, AUTOSOMAL = {
		for (b in nonpiv.offs) {
			mm = aperm.default(TR.MATR[fa_haps, mo_haps, haps[[b]], drop = F], perm = c(3, 1, 2), resize = TRUE) * probs[[b]]
			likel = likel * .colSums(mm, length(haps[[b]]), fa_len * mo_len)
		}
		if(pivtype==0) return(sum(likel))
		
		res = switch(pivtype, .rowSums(likel, fa_len, mo_len), .colSums(likel, fa_len, mo_len),
		{ T = TR.MATR[fa_haps, mo_haps, haps[[piv]], drop = F]
		  arr = as.vector(T) * as.vector(likel)
		  dim(arr) = dim(T)
		  .colSums(arr, fa_len * mo_len, length(haps[[piv]])) * probs[[piv]]
		})
		haps[[piv]] = haps[[piv]][res != 0]
		probs[[piv]] = res[res != 0]
		return(dat = list(probs = probs, haps = haps))
	}, 
	X = {
		offlik = rep.int(1, length(mo_haps))
		for (moff in nonpiv.offs[SEX[nonpiv.offs] == 1]) {
			mm = t.default(TR.MATR[[1]][mo_haps, haps[[moff]], drop = F]) * probs[[moff]]
			offlik = offlik * .colSums(mm, length(haps[[moff]]), mo_len)
		}
		likel = t.default(t.default(likel) * offlik)
		for (foff in nonpiv.offs[SEX[nonpiv.offs] == 2]) {
			mm = aperm.default(TR.MATR[[2]][fa_haps, mo_haps, haps[[foff]], drop = F], c(3, 1, 2), TRUE) * probs[[foff]]
			likel = likel * .colSums(mm, length(haps[[foff]]), fa_len * mo_len)
		}
		if(pivtype==0) return(sum(likel))
		
		res = switch(pivtype, .rowSums(likel, fa_len, mo_len), .colSums(likel, fa_len, mo_len),
		{ switch(SEX[piv],
			{  T = TR.MATR[[1]][mo_haps, haps[[piv]], drop = F]
			   res = .colSums(T * .colSums(likel, fa_len, mo_len), mo_len, length(haps[[piv]])) * probs[[piv]]
			},
			{  T = TR.MATR[[2]][fa_haps, mo_haps, haps[[piv]], drop = F]
			   arr = as.vector(T) * as.vector(likel)
			   dim(arr) = dim(T)
			   res = .colSums(arr, fa_len * mo_len, length(haps[[piv]])) * probs[[piv]]
			})
		})
		haps[[piv]] <- haps[[piv]][res != 0]
		probs[[piv]] <- res[res != 0]
		return(dat = list(probs = probs, haps = haps))
	})
}


.likelihoodSNP = function (x, marker, theta = NULL, afreq = NULL, logbase = NULL, TR.MATR = NULL, singleNum.geno = NULL) {
  	 if (inherits(x,'singleton')) stop("This function is not applicable to singleton objects")
    if (x$hasLoops) stop("Unbroken loops in pedigree.")
    if (is.null(theta) && is.null(TR.MATR) && is.null(x$model)) stop("No model set.")
    nInd = x$nInd;    ped = x$pedigree;    chrom = x$model$chrom
    if (is.null(singleNum.geno)) {
        if (length(marker) == 1) marker = x$markerdata[[marker]]
        if (max(marker) > 2) stop("Marker has more than two alleles. You should use 'likelihood_multiallele' instead.")
        if (is.null(afreq)) afreq = attr(marker, "afreq")
        singleNum.geno = .diallel2geno(marker)
    }
    if (is.null(afreq)) stop("Allele frequencies missing.")
    if (length(afreq) == 1) afreq = c(afreq, 1 - afreq)
    if (is.null(TR.MATR)) TR.MATR = .TRmatr(theta, chrom)
    
	switch(chrom, AUTOSOMAL = {
        hap_list <- list(`??` = 1:10, AA = 1:3, BB = 8:10, AB = 4:7, `A?` = 1:7, `B?` = 4:10)[singleNum.geno + 1]
        a = afreq[1];   b = afreq[2]
        init_p = x$initial_probs
        init_p[, x$founders] = init_p[, x$founders] * rep(c(a^2, 2*a*b, b^2), c(3, 4, 3))
        prob_list <- lapply(seq_len(nInd), function(i) init_p[, i][hap_list[[i]]])
    }, X = {
        haplo.poss_X <- list(list(`?` = 1:4, A = 1:2, B = 3:4, NULL, NULL, NULL), list(`??` = 1:10, AA = 1:3, BB = 8:10, AB = 4:7, `A?` = 1:7, `B?` = 4:10))
        hap_list <- lapply(seq_len(nInd), function(i) haplo.poss_X[[ped[i, "SEX"]]][[singleNum.geno[i] + 1]])
        a = afreq[1];  b = afreq[2]
        AfreqX <- list(male = c(a, a, b, b), female = rep(c(a^2, 2*a*b, b^2), c(3, 4, 3)))
        init_p_list = x$initial_probs
        for (i in x$founders) 
			init_p_list[[i]] <- init_p_list[[i]] * AfreqX[[ped[i, "SEX"]]]
        prob_list <- lapply(seq_len(nInd), function(i) init_p_list[[i]][hap_list[[i]]])
    })
    hap_list <- lapply(seq_len(nInd), function(i) hap_list[[i]][prob_list[[i]] != 0])
    prob_list <- lapply(prob_list, function(v) v[v != 0])
    if (is.null(dups <- x$loop_breakers)) {
        dat = list(probs = prob_list, haps = hap_list)
        for (sub in x$subnucs) dat = .peelSNP(dat, sub, SEX = ped[, "SEX"], chrom = chrom, TR.MATR)
        likelihood = dat
    }
    else {
        origs = match(dups[, 1], x$orig.ids)
        copies = match(dups[, 2], x$orig.ids)
        sumover = fast.grid(lapply(seq_along(origs), function(i) intersect(hap_list[[origs[i]]], hap_list[[copies[i]]])))
        sumoverlist = lapply(seq_len(nrow(sumover)), function(ri) sumover[ri, ])
        likelihood = 0
        for (r in sumoverlist) {
            probs = prob_list;  haps = hap_list
            for (i in seq_along(origs)) {
                haps[[origs[i]]] <- haps[[copies[i]]] <- r[i]
                probs[[origs[i]]] = probs[[origs[i]]][r[i] == hap_list[[origs[i]]]]
                probs[copies[i]] = list(1)
            }
            dat = list(probs = probs, haps = haps)
            for (sub in x$subnucs) dat = .peelSNP(dat, sub, SEX = ped[, "SEX"], chrom = chrom, TR.MATR)
            likelihood = likelihood + dat
        }
    }
    if (is.numeric(logbase)) log(likelihood, logbase) else likelihood
}


.TRmatr=function(theta, chrom) {
	if (is.null(theta)) stop("Argument 'theta' cannot be NULL.")
	haplo.single <- c('AD','AN','BD','BN')
	haplo.allpairs <- c('AADD','AADN','AANN','ABDD','ABDN','ABND','ABNN','BBDD','BBDN','BBNN')
	h <- c( c(1,.5,0,.5,.5*(1-theta),.5*theta,0,0,0,0), c(0,.5,1,0,.5*theta,.5*(1-theta),.5,0,0,0), c(0,0,0,.5,.5*theta,.5*(1-theta),0,1,.5,0), 
			c(0,0,0,0,.5*(1-theta),.5*theta,.5,0,.5,1) )
	dim(h) <- c(10,4); dimnames(h) <- list(haplo.allpairs, haplo.single)
	switch(chrom, 
	AUTOSOMAL = {
		T <- numeric(1000); dim(T) <- c(10,10,10); dimnames(T) <- list(haplo.allpairs, haplo.allpairs, haplo.allpairs)
		T[,,'AADD'] <- h[,'AD'] %*% t(h[,'AD'])
		T[,,'AADN'] <- h[,'AD'] %*% t(h[,'AN']) + h[,'AN'] %*% t(h[,'AD'])
		T[,,'AANN'] <- h[,'AN'] %*% t(h[,'AN'])
		T[,,'ABDD'] <- h[,'AD'] %*% t(h[,'BD']) + h[,'BD'] %*% t(h[,'AD'])
		T[,,'ABDN'] <- h[,'AD'] %*% t(h[,'BN']) + h[,'BN'] %*% t(h[,'AD'])
		T[,,'ABND'] <- h[,'AN'] %*% t(h[,'BD']) + h[,'BD'] %*% t(h[,'AN'])
		T[,,'ABNN'] <- h[,'AN'] %*% t(h[,'BN']) + h[,'BN'] %*% t(h[,'AN'])
		T[,,'BBDD'] <- h[,'BD'] %*% t(h[,'BD'])
		T[,,'BBDN'] <- h[,'BD'] %*% t(h[,'BN']) + h[,'BN'] %*% t(h[,'BD'])
		T[,,'BBNN'] <- h[,'BN'] %*% t(h[,'BN'])
		return(T)
	}, X = {
		TR_f <- numeric(400); dim(TR_f) <- c(4,10,10); dimnames(TR_f) <- list(haplo.single, haplo.allpairs, haplo.allpairs)
		TR_f['AD', , c('AADD','AADN','ABDD','ABDN')] <- h
		TR_f['AN', , c('AADN','AANN','ABND','ABNN')] <- h 
		TR_f['BD', , c('ABDD','ABND','BBDD','BBDN')] <- h 
		TR_f['BN', , c('ABDN','ABNN','BBDN','BBNN')] <- h
		return(list(h, TR_f))
	})
}


.diallel2geno <- function(marker) { #marker a numerical nInd * 2 matrix   
	#Coding genotypes as single integer: 00 -> 0, 11 -> 1, 22 -> 2, 12/21 -> 3, 01/10 -> 4, 02/20 -> 5.
	#Each pair of alleles is seen as an integer written in base 3, and this integer is permuted to fit with the above code.
	c(0,4,5,4,1,3,5,3,2)[colSums(c(3,1) * t(marker)) + 1]
}

.geno2diallel <- function(codedgenos) { #input: matrix of single-numerical genotypes. 
	#Ouput: Matrix with twice the number of columns, decoded as 1 -> 1 1, 2 -> 1 2, 3 -> 2 2, 4 -> 1 0, 5 -> 2 0
	decode = sapply(t(codedgenos+1),function(i) switch(i, c(0,0), c(1,1), c(2,2), c(1,2), c(1,0), c(2,0)))
	dim(decode) = c(2*ncol(codedgenos), nrow(codedgenos))
	t(decode)
}

