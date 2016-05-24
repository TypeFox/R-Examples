linkdat <-
function(ped, model=NULL, map=NULL, dat=NULL, freq=NULL, annotations=NULL, missing=0, header=FALSE, verbose=TRUE, ...) {  

	subnucs <- function(ped) {  	# output: peeling order of nuclear subfamilies. Format for each nuc: c(pivot,father,mother,offsp1,..), where pivot=0 for the last nuc.
		if (nrow(ped)==1) return(list())
      parents = unique(ped[, 2:3])
		parents = parents[ -match(0, parents[,1]), , drop=FALSE]
		list1 = lapply(nrow(parents):1, function(i) { 
			par = parents[i,]; 
			list(father = par[[1]], mother = par[[2]], offspring = as.vector(ped[,1])[ which(ped[,2]==par[[1]] & ped[,3]==par[[2]], useNames=FALSE) ]) } )  #listing all nucs
		res = list(); i = 1; k = 1
		while(length(list1) > 1)	{
			if (i > length(list1)) 	return(FALSE)
			sub = list1[[i]]
			subvec = unlist(sub)
			links = subvec[subvec %in% unlist(list1[-i])]
			if (length(links)==1) 	{
				res[[k]] <- c(sub, list(pivot = as.numeric(links), pivtype = match(links, c(sub[['father']], sub[['mother']]), nomatch=3)))
				list1 <- list1[-i]
				k <- k+1;	i <- 1
			} else i <- i+1
		}	
		res[[k]] <- c(list1[[1]], list(pivot = 0, pivtype = 0)) #final nuclear
		res
	}
	
    #if(!is.null(annotations) && (!is.null(map) || !is.null(dat))) stop("If 'annotations' is given, 'map' and 'dat' must be NULL")
	if(inherits(ped, c('linkdat', 'singleton'))) return(ped)	
	if(is.character(ped) && length(ped)==1) {
        skip = as.integer(header)
        first = scan(ped, what="", skip=skip, nlines=1, quiet=TRUE, ...)
        ncols = length(first)
        .numerical = !any(is.na(suppressWarnings(as.numeric(first))))
        ped = scan(ped, what="", skip=skip, quiet=TRUE, ...)
        ped = matrix(ped, ncol=ncols, byrow=TRUE)
    }
    ped = as.matrix(ped) # if ped is a data frame
    if(!exists(".numerical")) .numerical = !is.numeric(ped) # looks wrong, but makes sense in next line: If ped is a character matrix
    numerical = .numerical && !any(is.na(suppressWarnings(ped_num <- as.numeric(ped))))
    if(numerical) {
        dim(ped_num) = dim(ped)
        ped = ped_num
    }
   
    nrows = nrow(ped)
    if(nrows == 0) stop("Empty pedigree.")
   
    if(!is.null(map)) {
        if(is.null(dat)) stop("The 'map' and 'dat' arguments must either both be NULL, or both non-NULL")
        annotations = .readMap(map, dat, freq, verbose, numerical)
    }
    else if (is.character(annotations))
        annotations = list(alleles=switch(annotations, snp12=1:2, snpAB=c('A','B')), afreq=c(0.5,0.5))
      
    if (length(unique(ped[, 1])) < nrows | all(gsub(" ","", ped[, 2], fixed=TRUE) != 0)) { #added second test, in case all rows are singletons
        famids = unique.default(ped[, 1]) # these are characters here
        if(length(famids) > 1) 
			return(lapply(famids[order(as.numeric(famids))], function(fam) 
            linkdat(ped[ped[, 1]==fam, , drop=F], model=model, annotations=annotations, verbose=verbose, missing=missing)))
		else {
			famid = as.numeric(ped[1, 1])
			ped = ped[, -1, drop=F]
        }
	}
	else if (nrows == 1 && all(ped[, 3:4]==0)) { #singleton with famid!
        famid = as.numeric(ped[1, 1])
		ped = ped[, -1, drop=F]
    }
    else famid=1
	if (verbose) cat("Family ID: ", famid, ".\n", sep="")
	
    if(ncol(ped) < 5) stop("Too few columns: ID, FID, MID, SEX and AFF are mandatory.")
    pedcols = ped[, 1:5]
    if(!numerical) {
        if(any(grepl("[^0-9 ]", pedcols))) stop("Pedigree columns must be numeric.")
        pedcols = as.numeric(pedcols)
    }
    pedigree = matrix(pedcols, ncol=5, dimnames = list(NULL, c('ID', 'FID', 'MID', 'SEX', 'AFF')))
    .checkped(pedigree)
	orig.ids = as.vector(pedigree[, 1])
	nInd = nrows
	pedigree = relabel(pedigree, new=1:nInd)
	if (verbose) {
        if(nInd == 1) cat("Singleton ", ifelse(pedigree[, 'SEX'] == 1, 'male', 'female'), ", individual ID = ", orig.ids, '.\n', sep="")
        else cat(nInd, "individuals.\n")
        affs = sum(pedigree[, 'AFF'] == 2)
        if(affs > 0) cat(affs, "affected", nInd-affs, "non-affected.\n")
	}
	founders = as.integer(which(pedigree[,'FID'] == 0))
	nonfounders = as.integer(which(pedigree[,'FID'] > 0))
		
   #---peeling order of nuclear subfamilies---
	if (nInd>1) {
        subnucs = subnucs(pedigree)
        hasLoops = is.logical(subnucs) && !subnucs
        if(verbose) 
            if(hasLoops) cat("Loop(s) detected.\n")  
            else if(is.list(subnucs)) 
                cat(ant <- length(subnucs), "nuclear", ifelse(ant==1,"subfamily.\n", "subfamilies.\n"))
        class = 'linkdat'
    } else {
        class = c('singleton', 'linkdat')
        hasLoops = FALSE
        subnucs = NULL
    }
   
	#---creation of linkdat object---
	obj = structure(list(pedigree=pedigree, famid=famid, orig.ids=orig.ids, nInd=nInd, founders=founders, 
        nonfounders=nonfounders, hasLoops=hasLoops, subnucs=subnucs), class=class)
   
	#---adding markers---
	obj = setMarkers(obj, ped[, -(1:5), drop=F], annotations=annotations, missing=missing)
	if (verbose) if(obj$nMark == 1) cat("1 marker.\n") else cat(obj$nMark, "markers.\n")

	#----adding model----
 	if (!is.null(model)) obj = setModel(obj, model=model)

	if (verbose) cat('\n')
	invisible(obj)
}


singleton = function(id, sex=1, famid=1, verbose=FALSE, ...)  
    linkdat(ped = rbind(c(famid, id, 0, 0, sex, 1)), verbose=verbose, ...)

#print.singleton = function(x, ...) {
#   n = if(is.null(x$markerdata)) 0 else x$nMark 
#   cat("Singleton ", ifelse(x$pedigree[,'SEX']==1, 'male', 'female'), "\nID label: ", x$orig.ids, "\n", n, " markers\n", sep="")
#}

.checkped <- function(p) { #p a numeric matrix with 5 columns
    ID = p[,'ID']; FID = p[,'FID']; MID = p[,'MID']; SEX = p[,'SEX']; AFF=p[,'AFF']
   
   #singletons:
	if (nrow(p) == 1) {
        if(FID != 0 || MID != 0) stop("Singleton error: FID and MID must both be zero.")
        if(!SEX %in% 1:2) stop("Singleton error: Unknown sex.")
        return()
    }
    #real pedigrees:
	if (all(c(FID, MID)==0)) stop("Pedigree is not connected.")
    
    fatherErr = !FID %in% c(0,ID)
    motherErr = !MID %in% c(0,ID)
    self_ancest = sapply(seq_along(ID), function(i) ID[i] %in% ancestors(p, ID[i]))
	quick.check <- (all(SEX %in% 1:2) && all(AFF %in% 0:2) && all((FID>0)==(MID>0)) &&
                    !any(duplicated(ID)) && !any(fatherErr) && !any(motherErr ) && 
                    all(SEX[match(FID[FID!=0], ID)] == 1) && 
                    all(SEX[match(MID[MID!=0], ID)] == 2) && 
					!any(self_ancest))
    
    if (quick.check)  #if all tests are passed
		return()
	else {
        for (i in seq_along(ID)) {
            if (!SEX[i] %in% 1:2) cat("Individual ", ID[i],": SEX must be either 1 (male) or 2 (female).\n", sep="")
            if (!AFF[i] %in% 0:2) cat("Individual ", ID[i],": Affection status must be either 0 (unknown), 1 (non-affected) or 2 (affected).\n", sep="")
            if ((FID[i]>0) != (MID[i]>0)) cat("Individual ", ID[i],": Only one parent in the pedigree is not allowed. Either both parents or none must be specified.\n", sep="")
            if (i>1 && ID[i] %in% ID[1:(i-1)]) {cat("Individual ", ID[i],": ID not unique.\n", sep=""); next}
            if (fatherErr[i]) cat("Individual ", ID[i],": Father's ID (",FID[i],") does not appear in ID column.\n", sep="")
            else if (FID[i] != 0 && SEX[match(FID[i], ID)] != 1) cat("Individual ", ID[i],": Father is not male.\n", sep="")
            if (motherErr[i]) cat("Individual ", ID[i],": Mother's ID (",MID[i],") does not appear in ID column.\n", sep="")
            else if (MID[i] != 0 && SEX[match(MID[i], ID)] != 2) cat("Individual ", ID[i],": Mother is not female.\n", sep="")
            if (self_ancest[i]) cat("Individual ", ID[i]," is ", switch(SEX[i],"his", "her"), " own ancestor.\n", sep="")
        }
        stop("Pedigree errors detected.")
	}
}
