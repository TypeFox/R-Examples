
################################
### Extending to 3 and 4 alleles
################################

.peelFAST <- function(dat, sub, SEX, chrom, TR.MATR) {
	probs = dat[["probs"]];   haps = dat[["haps"]]   
   father = sub[["father"]]; mother = sub[["mother"]]
   fa_haps = haps[[father]];   mo_haps = haps[[mother]]
   fa_len = length(fa_haps);   mo_len = length(mo_haps);  famo_prod = fa_len*mo_len
	likel = probs[[father]] %*% t.default(probs[[mother]])
	piv = sub[["pivot"]]; pivtype = sub[["pivtype"]]
	offs = sub[['offspring']]; nonpiv.offs = offs[offs != piv]

	switch(chrom, 
   AUTOSOMAL = {
		for (b in nonpiv.offs) {
         b_haps = haps[[b]]; b_probs = probs[[b]]
         mm = TR.MATR[fa_haps, mo_haps, b_haps, drop = F]
         if(prod(b_probs) < 1) 
            mm = mm * rep(b_probs, each=famo_prod)
			likel = likel * .rowSums(mm, famo_prod, length(b_haps))
		}
		if(pivtype==0) return(sum(likel))
		
      piv_haps = haps[[piv]]
		res = switch(pivtype, .rowSums(likel, fa_len, mo_len), .colSums(likel, fa_len, mo_len),
		{ 
        T = TR.MATR[fa_haps, mo_haps, piv_haps, drop = F]
		  arr = as.numeric(T) * as.numeric(likel)
		  .colSums(arr, famo_prod, length(piv_haps)) * probs[[piv]]
		})
		dat$haps[[piv]] = piv_haps[res != 0]
		dat$probs[[piv]] = res[res != 0]
      if(sum(res)==0) attr(dat, 'impossible') = TRUE
		return(dat)
	}, 
	X = {
		offlik = rep.int(1, length(mo_haps))
		for (moff in nonpiv.offs[SEX[nonpiv.offs] == 1]) {
			mm = t.default(TR.MATR[[1]][mo_haps, haps[[moff]], drop = F]) * probs[[moff]]
			offlik = offlik * .colSums(mm, length(haps[[moff]]), mo_len)
		}
		likel = t.default(t.default(likel) * offlik)
		
      for (foff in nonpiv.offs[SEX[nonpiv.offs] == 2]) {
			foff_haps = haps[[foff]]; foff_probs = probs[[foff]]
         mm = TR.MATR[[2]][fa_haps, mo_haps, foff_haps, drop = F]
         if(prod(foff_probs) < 1) 
            mm = mm * rep(foff_probs, each=famo_prod)
			likel = likel * .rowSums(mm, famo_prod, length(foff_haps))
         
         #mm = aperm.default(TR.MATR[[2]][fa_haps, mo_haps, haps[[foff]], drop = F], c(3, 1, 2), TRUE) * probs[[foff]]
			#likel = likel * .colSums(mm, length(haps[[foff]]), fa_len * mo_len)
		}
		if(pivtype==0) return(sum(likel))
		
      piv_haps = haps[[piv]]
		res = switch(pivtype, .rowSums(likel, fa_len, mo_len), .colSums(likel, fa_len, mo_len),
         {switch(SEX[piv],
            { T = TR.MATR[[1]][mo_haps, piv_haps, drop = F]
			     .colSums(T * .colSums(likel, fa_len, mo_len), mo_len, length(piv_haps)) * probs[[piv]]
            },
            { T = TR.MATR[[2]][fa_haps, mo_haps, piv_haps, drop = F]
			     arr = as.numeric(T) * as.numeric(likel)
			     .colSums(arr, famo_prod, length(piv_haps)) * probs[[piv]]
            }
         )})
         
		dat$haps[[piv]] = piv_haps[res > 0]
		dat$probs[[piv]] = res[res > 0]
      if(sum(res)==0) attr(dat, 'impossible') = TRUE
		return(dat)
	})
}


.initialCalc = function(x, afreq, chrom) {
   d = x$model$dfreq
   N = length(afreq)
   G = N*(N+1)/2 # number of marker genotypes
   Ghet = N*(N-1)/2  # number of heterozygous genotypes
   gtMatr = allGenotypes(N)  # matrix with 2 columns and G rows
	
   gtnames = paste(letters[gtMatr[,1]], letters[gtMatr[,2]], sep="")
   haplonames = paste(rep(gtnames[1:N], each=3), c('DD','DN','NN'), sep="")
   if (N > 1) 
      haplonames = c(haplonames, paste(rep(gtnames[-(1:N)], each=4), c('DD','DN','ND','NN'), sep=""))
   haps_init <- list(`??` = seq_along(haplonames))
   for(gt in gtnames) 
      haps_init[[gt]] = grep(gt, haplonames)
   for(a in letters[1:N]) 
      haps_init[[paste(a,'?',sep="")]] = grep(a, haplonames)
   
   penets = x$model$penetrances
   if(chrom=="X") penets = penets$female
   p = penets[c(rep.int(3:1, N), rep.int(c(3,2,2,1), Ghet))]  # P(aff | geno). Note that P(non-aff | geno) = 1-p
   penmatr = cbind(1, 1-p, p) #rownames = list(c('AADD','AADN','AANN','BB',... 'AB', ... 'AC',...
	
   disFreq = c(rep(c(d^2, 2*d*(1-d), (1-d)^2), N), rep(c(d^2, d*(1-d), d*(1-d), (1-d)^2), Ghet))
   
   gtFreq = afreq[gtMatr[,1]] * afreq[gtMatr[,2]] * rep.int(1:2, c(N, Ghet))
   gtFreqExtended = rep(gtFreq, rep.int(3:4, c(N, Ghet)))
   
   aff = x$pedigree[, 'AFF']
   switch(chrom,
	AUTOSOMAL = {
		initial_probs = penmatr[, aff + 1]
		initial_probs[, x$founders] = initial_probs[, x$founders] * disFreq * gtFreqExtended
      #dimnames(initial_probs) = list(haplonames, x$orig.ids) # kan sloyfes
   },
   X = {
   haps_init_male = list(`?` = seq_len(2*N))
   haps_init_male[seq_len(N)+1] = lapply(seq_len(N)*2, function(i) c(i-1, i))
   haps_init = list(male = haps_init_male, female=haps_init)
   
   pM = x$model$penetrances$male[rep.int(2:1, N)] #P(aff | geno) for males
   penmatrM = cbind(1, 1-pM, pM) # dimnames = list(c('AD','AN','BD','BN', ...), 1:3))
	
   penmatrX = list(male=penmatrM, female=penmatr)
	disFreqX = list(male=rep(c(d, 1-d), N), female=disFreq)
   gtFreqX = list(male=rep(afreq, each=2), female=gtFreqExtended)
   
   sex = x$pedigree[, 'SEX']
   initial_probs = lapply(1:x$nInd, function(i) penmatrX[[ sex[i] ]][, aff[i] + 1])
	for (i in x$founders) 	
      initial_probs[[i]] = initial_probs[[i]] * disFreqX[[ sex[i] ]] * gtFreqX[[ sex[i] ]]
   })
   
   list(initial_probs=initial_probs, haps_init=haps_init)
}
    
likelihood_LINKAGE = function (x, marker, theta = NULL, afreq = NULL, logbase = NULL, 
                           TR.MATR = NULL, initialCalc = NULL, singleNum.geno = NULL) {
  	 if (inherits(x,'singleton')) stop("This function is not applicable to singleton objects")
    if (x$hasLoops) stop("Unbroken loops in pedigree.")
    if (is.null(x$model)) stop("No model set.")
    
    nInd = x$nInd
    ped = x$pedigree
    chrom = x$model$chrom
    SEX = ped[, "SEX"]
    if (is.null(singleNum.geno)) {
        if (length(marker) == 1) marker = x$markerdata[[marker]]
        if (max(marker) > 4) stop("Marker has more than 4 alleles. You should use 'likelihood()' instead.")
        if (is.null(afreq)) afreq = attr(marker, "afreq")
        singleNum.geno = .diallel2genoNEW(marker, n=length(afreq))
    }
    
    if (is.null(TR.MATR)) TR.MATR = .TRmatrNEW(theta, length(afreq), chrom)
    if (is.null(initialCalc)) initialCalc = .initialCalc(x, afreq, chrom)
    
    init_probs = initialCalc$initial_probs
    
    switch(chrom, AUTOSOMAL = {
        hap_list = initialCalc$haps_init[singleNum.geno + 1]
        prob_list = lapply(seq_len(nInd), function(i) init_probs[hap_list[[i]], i])
    }, X = {
        haps_init = initialCalc$haps_init
        hap_list = lapply(seq_len(nInd), function(i) haps_init[[SEX[i]]][[singleNum.geno[i] + 1]])
        prob_list = lapply(seq_len(nInd), function(i) init_probs[[i]][hap_list[[i]]])
    })
    hap_list = lapply(seq_len(nInd), function(i) hap_list[[i]][prob_list[[i]] > 0])
    prob_list = lapply(prob_list, function(v) v[v > 0])
    dat = list(probs = prob_list, haps = hap_list)
    attr(dat, 'impossible') = FALSE
    
    if (is.null(dups <- x$loop_breakers)) {
      for (sub in x$subnucs) {
         dat = .peelFAST(dat, sub, SEX = SEX, chrom = chrom, TR.MATR)
         if (sub$pivtype > 0 && attr(dat, 'impossible')) return(ifelse(is.numeric(logbase), -Inf, 0))
      }
      likelihood = dat
    } 
    else {
        origs = match(dups[, 1], x$orig.ids)
        copies = match(dups[, 2], x$orig.ids)
        loopgrid = fast.grid(lapply(seq_along(origs), function(i)  { 
                  seq_along(hap_list[[origs[i]]])[hap_list[[origs[i]]] %in% hap_list[[copies[i]]]] 
        }), as.list=TRUE)
        likelihood = 0
        for (r in loopgrid) {
            haps = dat$haps; probs = dat$probs
            for (i in seq_along(origs)) {
                orig = origs[i]; copy = copies[i]
                haps[[orig]] <- haps[[copy]] <- haps[[orig]][r[i]]
                probs[[orig]] = probs[[orig]][r[i]]
                if(sum(probs[[orig]])==0) print("Loop-loekke: Alle sannsynligheter er null. Magnus lurer paa om dette kan gi feilmelding.")
                probs[copy] = list(1)
            }
            dat1 = list(probs = probs, haps = haps); attr(dat1, 'impossible') = FALSE
            for (sub in x$subnucs) {
               dat1 = .peelFAST(dat1, sub, SEX = SEX, chrom = chrom, TR.MATR)
               if (sub$pivtype > 0 && attr(dat1, 'impossible')) {break} #if impossible data - break out of ES-algorithm and go to next r in loopgrid.
               if (sub$pivtype == 0) likelihood = likelihood + dat1
            }
            
        }
    }
    if (is.numeric(logbase)) log(likelihood, logbase) else likelihood
}

.haplonames = function(n) { # TODO Unngaa paste.....tar mye tid i lod().
   gtMatr = allGenotypes(n)
   gt = paste0(letters[gtMatr[,1]], letters[gtMatr[,2]])
   dishom = c('DD','DN','NN'); dishet = c('DD','DN','ND','NN')
   res = paste0(rep(gt[1:n], each=3), dishom)
   if(n > 1) 
      res = c(res, paste0(rep(gt[-(1:n)], each=4), dishet))
   res
}

.TRmatrNEW = function(theta, n, chrom=c('AUTOSOMAL', 'X')) { 
   if(!is.numeric(theta)) cat('theta',theta,'hva var det?')
   stopifnot(is.numeric(theta), length(theta)==1, theta >= 0, is.numeric(n), length(n)==1, n > 0)
   # n = #marker alleles
   k = n*(2*n + 1) # = number of haplotype pairs with marker and disease locus.
   haplo.single = paste0(rep(letters[1:n], each=2), c('D','N'))
   haplo.allpairs = .haplonames(n) 
   #n=3 --> c('AADD','AADN','AANN','BBDD','BBDN','BBNN', 'CCDD','CCDN','CCNN', 'ABDD','ABDN','ABND','ABNN', 'ACDD','ACDN','ACND','ACNN', 'BCDD','BCDN','BCND','BCNN') 
   homD = c(1, .5, 0); homN = c(0, .5, 1); hom0 = c(0,0,0)
   het1D = c(.5, .5*(1-theta), .5*theta, 0)# e.g. ABDD, ABDN, ABND, ABNN -> AD 
   het2D = c(.5, .5*theta, .5*(1-theta), 0)# e.g. ABDD, ABDN, ABND, ABNN -> BD 
   het1N = c(0, .5*theta, .5*(1-theta), 0.5)# e.g. ABDD, ABDN, ABND, ABNN -> AN
   het2N = c(0, .5*(1-theta), .5*theta, 0.5)# e.g. ABDD, ABDN, ABND, ABNN -> BN 
   het0 = c(0,0,0,0)
   
   ### computing the k * 2n matrix h, containing haplotype transmission probabilities from 1 parent. 
   if(n == 1)
      h = cbind( 
         aD = c(homD),
         aN = c(homN))
   else if(n == 2)
      h = cbind( 
         aD = c(homD, hom0, het1D),
         aN = c(homN, hom0, het1N),
         bD = c(hom0, homD, het2D),
         bN = c(hom0, homN, het2N))
   else if(n == 3)
      h = cbind( 
         aD = c(homD, hom0, hom0, het1D, het1D, het0),
         aN = c(homN, hom0, hom0, het1N, het1N, het0),
         bD = c(hom0, homD, hom0, het2D, het0, het1D),
         bN = c(hom0, homN, hom0, het2N, het0, het1N),
         cD = c(hom0, hom0, homD, het0, het2D, het2D),
         cN = c(hom0, hom0, homN, het0, het2N, het2N))
   else if(n == 4)
      h = cbind( 
         aD = c(homD, hom0, hom0, hom0, het1D, het1D, het1D, het0, het0, het0),
         aN = c(homN, hom0, hom0, hom0, het1N, het1N, het1N, het0, het0, het0),
         bD = c(hom0, homD, hom0, hom0, het2D, het0, het0, het1D, het1D, het0),
         bN = c(hom0, homN, hom0, hom0, het2N, het0, het0, het1N, het1N, het0),
         cD = c(hom0, hom0, homD, hom0, het0, het2D, het0, het2D, het0, het1D),
         cN = c(hom0, hom0, homN, hom0, het0, het2N, het0, het2N, het0, het1N),
         dD = c(hom0, hom0, hom0, homD, het0, het0, het2D, het0, het2D, het2D),
         dN = c(hom0, hom0, hom0, homN, het0, het0, het2N, het0, het2N, het2N))  
   else stop("More than 4 alleles not implemented yet.")
   
   chrom = match.arg(chrom)
   switch(chrom, 
   AUTOSOMAL = {
      T <- numeric(k^3); dim(T) <- rep(k, 3); dimnames(T) <- list(haplo.allpairs, haplo.allpairs, haplo.allpairs)
      T[,,'aaDD'] <- h[,'aD'] %*% t.default(h[,'aD'])
      T[,,'aaDN'] <- h[,'aD'] %*% t.default(h[,'aN']) + h[,'aN'] %*% t.default(h[,'aD'])
      T[,,'aaNN'] <- h[,'aN'] %*% t.default(h[,'aN'])
      if(n >= 2) {
         T[,,'bbDD'] <- h[,'bD'] %*% t.default(h[,'bD'])
         T[,,'bbDN'] <- h[,'bD'] %*% t.default(h[,'bN']) + h[,'bN'] %*% t.default(h[,'bD'])
         T[,,'bbNN'] <- h[,'bN'] %*% t.default(h[,'bN'])
         T[,,'abDD'] <- h[,'aD'] %*% t.default(h[,'bD']) + h[,'bD'] %*% t.default(h[,'aD'])
         T[,,'abDN'] <- h[,'aD'] %*% t.default(h[,'bN']) + h[,'bN'] %*% t.default(h[,'aD'])
         T[,,'abND'] <- h[,'aN'] %*% t.default(h[,'bD']) + h[,'bD'] %*% t.default(h[,'aN'])
         T[,,'abNN'] <- h[,'aN'] %*% t.default(h[,'bN']) + h[,'bN'] %*% t.default(h[,'aN'])
      }
      if(n >= 3) {
         T[,,'ccDD'] <- h[,'cD'] %*% t.default(h[,'cD'])
         T[,,'ccDN'] <- h[,'cD'] %*% t.default(h[,'cN']) + h[,'cN'] %*% t.default(h[,'cD'])
         T[,,'ccNN'] <- h[,'cN'] %*% t.default(h[,'cN'])
         T[,,'acDD'] <- h[,'aD'] %*% t.default(h[,'cD']) + h[,'cD'] %*% t.default(h[,'aD'])
         T[,,'acDN'] <- h[,'aD'] %*% t.default(h[,'cN']) + h[,'cN'] %*% t.default(h[,'aD'])
         T[,,'acND'] <- h[,'aN'] %*% t.default(h[,'cD']) + h[,'cD'] %*% t.default(h[,'aN'])
         T[,,'acNN'] <- h[,'aN'] %*% t.default(h[,'cN']) + h[,'cN'] %*% t.default(h[,'aN'])
         T[,,'bcDD'] <- h[,'bD'] %*% t.default(h[,'cD']) + h[,'cD'] %*% t.default(h[,'bD'])
         T[,,'bcDN'] <- h[,'bD'] %*% t.default(h[,'cN']) + h[,'cN'] %*% t.default(h[,'bD'])
         T[,,'bcND'] <- h[,'bN'] %*% t.default(h[,'cD']) + h[,'cD'] %*% t.default(h[,'bN'])
         T[,,'bcNN'] <- h[,'bN'] %*% t.default(h[,'cN']) + h[,'cN'] %*% t.default(h[,'bN'])
      }
      if(n >= 4) {
         T[,,'ddDD'] <- h[,'dD'] %*% t.default(h[,'dD'])
         T[,,'ddDN'] <- h[,'dD'] %*% t.default(h[,'dN']) + h[,'dN'] %*% t.default(h[,'dD'])
         T[,,'ddNN'] <- h[,'dN'] %*% t.default(h[,'dN'])
         T[,,'adDD'] <- h[,'aD'] %*% t.default(h[,'dD']) + h[,'dD'] %*% t.default(h[,'aD'])
         T[,,'adDN'] <- h[,'aD'] %*% t.default(h[,'dN']) + h[,'dN'] %*% t.default(h[,'aD'])
         T[,,'adND'] <- h[,'aN'] %*% t.default(h[,'dD']) + h[,'dD'] %*% t.default(h[,'aN'])
         T[,,'adNN'] <- h[,'aN'] %*% t.default(h[,'dN']) + h[,'dN'] %*% t.default(h[,'aN'])
         T[,,'bdDD'] <- h[,'bD'] %*% t.default(h[,'dD']) + h[,'dD'] %*% t.default(h[,'bD'])
         T[,,'bdDN'] <- h[,'bD'] %*% t.default(h[,'dN']) + h[,'dN'] %*% t.default(h[,'bD'])
         T[,,'bdND'] <- h[,'bN'] %*% t.default(h[,'dD']) + h[,'dD'] %*% t.default(h[,'bN'])
         T[,,'bdNN'] <- h[,'bN'] %*% t.default(h[,'dN']) + h[,'dN'] %*% t.default(h[,'bN'])
         T[,,'cdDD'] <- h[,'cD'] %*% t.default(h[,'dD']) + h[,'dD'] %*% t.default(h[,'cD'])
         T[,,'cdDN'] <- h[,'cD'] %*% t.default(h[,'dN']) + h[,'dN'] %*% t.default(h[,'cD'])
         T[,,'cdND'] <- h[,'cN'] %*% t.default(h[,'dD']) + h[,'dD'] %*% t.default(h[,'cN'])
         T[,,'cdNN'] <- h[,'cN'] %*% t.default(h[,'dN']) + h[,'dN'] %*% t.default(h[,'cN'])
      }
      return(T)
   }, 
   X = { 
      TR_f <- numeric(2*n*k^2); dim(TR_f) <- c(2*n,k,k); dimnames(TR_f) <- list(haplo.single, haplo.allpairs, haplo.allpairs)
      if(n==1) {
         TR_f['aD', , c('aaDD', 'aaDN')] <- h
         TR_f['aN', , c('aaDN', 'aaNN')] <- h
      }
      else if(n == 2) {
         TR_f['aD', , c('aaDD','aaDN','abDD','abDN')] <- h
         TR_f['aN', , c('aaDN','aaNN','abND','abNN')] <- h 
         TR_f['bD', , c('abDD','abND','bbDD','bbDN')] <- h 
         TR_f['bN', , c('abDN','abNN','bbDN','bbNN')] <- h
      }
      else if(n == 3) {
         TR_f['aD', , c('aaDD','aaDN','abDD','abDN','acDD','acDN')] <- h 
         TR_f['aN', , c('aaDN','aaNN','abND','abNN','acND','acNN')] <- h 
         TR_f['bD', , c('abDD','abND','bbDD','bbDN','bcDD','bcDN')] <- h 
         TR_f['bN', , c('abDN','abNN','bbDN','bbNN','bcND','bcNN')] <- h 
         TR_f['cD', , c('acDD','acND','bcDD','bcND','ccDD','ccDN')] <- h 
         TR_f['cN', , c('acDN','acNN','bcDN','bcNN','ccDN','ccNN')] <- h 
      }
      else if(n == 4) {
         TR_f['aD', , c('aaDD','aaDN','abDD','abDN','acDD','acDN','adDD','adDN')] <- h 
         TR_f['aN', , c('aaDN','aaNN','abND','abNN','acND','acNN','adND','adNN')] <- h 
         TR_f['bD', , c('abDD','abND','bbDD','bbDN','bcDD','bcDN','bdDD','bdDN')] <- h 
         TR_f['bN', , c('abDN','abNN','bbDN','bbNN','bcND','bcNN','bdND','bdNN')] <- h 
         TR_f['cD', , c('acDD','acND','bcDD','bcND','ccDD','ccDN','cdDD','cdDN')] <- h 
         TR_f['cN', , c('acDN','acNN','bcDN','bcNN','ccDN','ccNN','cdND','cdNN')] <- h 
         TR_f['dD', , c('adDD','adND','bdDD','bdND','cdDD','cdND','ddDD','ddDN')] <- h 
         TR_f['dN', , c('adDN','adNN','bdDN','bdNN','cdDN','cdNN','ddDN','ddNN')] <- h 
      }
      return(list(male=h, female=TR_f))
   })
}



.diallel2genoNEW <- function(marker, n) { #marker: a numerical nInd * 2 matrix   
   # Coding genotypes as single integer by interpreting the allele pair as an integer written in base n+1.
   if (n == 1) {
      # Target code 00=0, 11=1, 01/10=2  
      # Base 2: 00=0, 01=1, 10=2, 11=3
      perm = c(0, 2, 2, 1)
   } else if(n == 2) { 
      # Target code 00=0, 11=1, 22=2, 12/21=3, 01/10=4, 02/20=5  
      # Base 3: 00=0, 01=1, 02=2, 10=3, 11=4, 12=5, 20=6, 21=7, 22=8
      perm = c(0, 4, 5, 4, 1, 3, 5, 3, 2)
   } else if(n == 3) {
      # Target: 00=0, 11=1, 22=2, 33=3, 12/21=4, 13/31=5, 23/32=6, 01/10=7, 02/20=8, 03/30=9  
      # Base 4: 00=0, 01=1, 02=2, 03=3, 10=4, 11=5, 12=6, 13=7, 20=8, 21=9, 22=10, 23=11, 30=12, 31=13, 32=14, 33=15
      perm = c(0,7,8,9,7,1,4,5,8,4,2,6,9,5,6,3)
   } else if(n == 4) {
      # Target: 00=0, 11=1, 22=2, 33=3, 44=4, 12/21=5, 13/31=6, 14/41=7, 23/32=8, 24/42=9, 34/43=10, 01/10=11, 02/20=12, 03/30=13, 04/40=14  
      # Base 5: 00=0, 01=1, 02=2, 03=3, 04=4, 10=5, 11=6, 12=7, 13=8, 14=9, 20=10, 21=11, 22=12, 23=13, 24=14, 30=15, 31=16, 32=17, 33=18, 34=19, 40=20, 41=21, 42=22, 43=23, 44=24
      perm = c(0,11,12,13,14,11,1,5,6,7,12,5,2,8,9,13,6,8,3,10,14,7,9,10,4)
   }
   perm[marker[, 1]*(n+1) + marker[, 2] + 1]
}

.geno2diallelNEW <- function(codedgenos, n) { #input: matrix of single-numerical genotypes. 
   codedgenos = as.matrix(codedgenos)
   ncols = ncol(codedgenos)
   #Ouput: Matrix with twice the number of columns, decoded as 1 -> 1 1, 2 -> 1 2, 3 -> 2 2, 4 -> 1 0, 5 -> 2 0
   if(n==1)
      gt=matrix(c(0,0, 1,1, 1,0), nrow=2)
   else if(n==2)
      gt=matrix(c(0,0, 1,1, 2,2, 1,2, 1,0, 2,0), nrow=2)
   else if(n==3)
      gt=matrix(c(0,0, 1,1, 2,2, 3,3, 1,2, 1,3, 2,3, 1,0, 2,0, 3,0), nrow=2)
   else if(n==4)
      gt=matrix(c(0,0, 1,1, 2,2, 3,3, 4,4, 1,2, 1,3, 1,4, 2,3, 2,4, 3,4, 1,0, 2,0, 3,0, 4,0), nrow=2)
   res = matrix(numeric(), nrow=nrow(codedgenos), ncol=2*ncols)
   even = 2*seq_len(ncols)
   alleles = as.vector(gt[, codedgenos+1])
   res[, even-1] = alleles[seq(1, length(alleles), by=2)]
   res[, even]   = alleles[seq(2, length(alleles), by=2)]
   res
}
