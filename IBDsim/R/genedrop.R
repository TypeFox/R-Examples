genedrop <-
function(x, map, condition=NULL, model="chi", skip.recomb=NULL) { #x=linkdat object
	ped = x$pedigree
    chrom = attr(map, "chromosome")
    #if(chrom==23) chrom="X"
    #if(!is.null(x$model$chrom) && chrom != x$model$chrom) stop(sprintf("Map chromosome = %s, but disease model chromosome = %s.", chrom, x$model$chrom))
    
    h = distribute.founder.alleles(x, chrom)
	nonfounders.order = dropping.order(x)
    
	if(is.null(condition)) {
		if(chrom < 23)
			for (i in nonfounders.order) {
                parents = ped[i, c('FID', 'MID')]
				h[[i]] = lapply(1:2, function(k) meiosis(h[[ parents[k] ]], map=map[[k]], model=model, skip.recomb=parents[k] %in% skip.recomb))
			}
		else if(chrom==23)
			for (i in nonfounders.order) {
				father = ped[i, 'FID']; mother = ped[i, 'MID']
				maternal.gamete = meiosis(h[[mother]], map=map$female, model=model, skip.recomb=mother %in% skip.recomb)
				if(ped[i,'SEX']==1) 	h[[i]] = list(maternal.gamete, maternal.gamete)
				else					h[[i]] = list(h[[father]][[1]], maternal.gamete)
			}	
		else stop("Chromosome is neither NULL, 'AUTOSOMAL' or 'X'.")
	}
	else {
		zero = condition$'0'; one = condition$'1'; two = condition$'2'; atm1 = condition$'atmost1';
		dis.fou = one[one %in% x$founders]; if(length(dis.fou)!=1) stop("Obligate carriers must include exactly 1 founder.")
		dis.al = h[[dis.fou]][[1]][[2]] #h[[dis.fou]][[1]] is matrix with 1 row.
		dis.locus = runif(1, min=0, max=attr(map, 'length_Mb'))
	
		carry_code = lapply(ped[,'ID'], function(j) match(T, c(j %in% zero, j %in% one, j %in% two, j %in% atm1), nomatch=5))
		COND = list(c(locus = dis.locus, allele = dis.al, action = 1), c(locus = dis.locus, allele = dis.al, action = 2), NULL) #action: 1=force; 2=avoid
		if(chrom < 23)
            for (i in nonfounders.order) {
				parents = ped[i, c('FID', 'MID')]; skip = parents %in% skip.recomb;	
				condits = COND[.decide_action(dis.al, dis.locus, h[parents], carry_code[[i]])] #returns list of 2 elements
				h[[i]] = lapply(1:2, function(k) 
					meiosis(h[[ parents[k] ]], map=map[[k]], model=model, skip.recomb=skip[k], condition = condits[[k]]))
		} 
		else stop("X-linked conditional genedrop is not implemented yet.")
		attr(h, 'dis.locus') = dis.locus; attr(h, 'dis.allele') = dis.al
	}
	attr(h, "length_Mb") <- attr(map, "length_Mb");	attr(h, "chromosome") <- chrom
	h
}

dropping.order = function(x) {
    # output: vector of all nonfounders, ordered such that children always come after their parents.
    taken = x$founders
    remaining = x$nonfounders
    while(length(remaining) > 0)
        for(k in seq_along(remaining)) {
            id = remaining[k]
            if(all(parents(x, id) %in% taken)) {
                taken = c(taken, id)
                remaining = remaining[-k]
                break
            }
        }
    setdiff(taken, x$founders)
}

distribute.founder.alleles = function(x, chrom="AUTOSOMAL") {
	ped = x$pedigree; fou = x$founders; h = vector("list", x$nInd)
	if(is.null(chrom)) chrom="AUTOSOMAL"
    if(is.numeric(chrom)) chrom=ifelse(chrom<23, "AUTOSOMAL", "X")
    
    if(chrom=="AUTOSOMAL")
		aux = cbind(rep.int(0, 2*length(fou)), seq_len(2*length(fou)))
	else {
		alleles = numeric(nfou <- length(fou)); k=1; 
		for(i in seq_along(fou)) {sex=ped[fou[i], 'SEX']; alleles[c(2*i-1,2*i)] = c(k, k + sex -1); k=k+sex}
		aux = cbind(rep.int(0, 2*nfou), alleles, deparse.level=0)
	}	
	h[fou] <- lapply(2*seq_along(fou), function(i) list(aux[i-1,,drop=F], aux[i,,drop=F]))
	h
}
	
.decide_action <- 
function(dis.al, dis.locus, parental.haplos, carry_code) { #carry_code = 1 (0 disease alleles), 2 (1), 3 (2) or 4 (at most 1);  
	fa = sum(dis.al == .getAlleles(parental.haplos[[1]], posvec=dis.locus)) 
	mo = sum(dis.al == .getAlleles(parental.haplos[[2]], posvec=dis.locus))
	f = 1; a = 2; nocond = 3; co = c(nocond,nocond)   #f = "force"; a = "avoid"
	err = function() stop("Impossible condition.")
	
    if(fa==0 && mo==0) switch(carry_code, NULL, err(), err(), NULL)
	else if(fa==1 && mo==0) switch(carry_code, co[1]<-a, co[1]<-f, err(), co[1]<-sample.int(2,1))
	else if(fa==0 && mo==1) switch(carry_code, co[2]<-a, co[2]<-f, err(), co[2]<-sample.int(2,1))
	else if(fa==2 && mo==0) switch(carry_code, err(), NULL, err(), NULL)
	else if(fa==0 && mo==2) switch(carry_code, err(), NULL, err(), NULL)
	else if(fa==1 && mo==1) switch(carry_code, co<-c(a,a), co<-sample.int(2), co<-c(f,f), co<-sample(list(c(f,a), c(a,f), c(a,a)), size=1)[[1]])
	else if(fa==1 && mo==2) switch(carry_code, err(), co[1]<-a, co[1]<-f, co[1]<-a)
	else if(fa==2 && mo==1) switch(carry_code, err(), co[2]<-a, co[2]<-f, co[2]<-a)
	else if(fa==2 && mo==2) switch(carry_code, err(), err(), co[2]<-f, co[2]<-a)
	return(co)
}

getLocus = function(x, h, locus) {
	marker = t.default(sapply(h, .getAlleles, posvec=locus))
	setMarkers(x, marker)
}

