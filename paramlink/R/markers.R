addMarker = function(x, m, alleles=NULL, afreq=NULL, missing=0) {
    if(is.matrix(m) || is.data.frame(m)) stopifnot(nrow(m)==x$nInd, ncol(m)==2)
    if(inherits(m, "marker"))
        m = list(m)
    if(is.list(m) && all(sapply(m, inherits, what="marker"))) 
        return(setMarkers(x, structure(c(x$markerdata, m), class="markerdata")))
    if(!is.list(m) && length(m)==1)    
        m = matrix(m, ncol=2, nrow=x$nInd)  #gives a nice way of setting an empty or everywhere-homozygous marker, e.g.: x=addMarker(x,0)
    mm = .createMarkerObject(m, alleles, afreq=afreq, missing=missing)
    setMarkers(x, structure(c(x$markerdata, list(mm)), class="markerdata"))
}

setMarkers <- function(x, m, annotations = NULL, missing=0) {
    if(is.null(m))                                markerdata_list = NULL
    else if(inherits(m, "marker"))                markerdata_list = structure(list(m), class="markerdata")
    else if(is.list(m) && all(sapply(m, inherits, what="marker"))) markerdata_list = structure(m, class="markerdata")
    else if(inherits(m, "markerdata"))            markerdata_list = m
    else if ((n <- ncol(m <- as.matrix(m))) == 0)    markerdata_list = NULL
    else {
        if(is.character(m[1,1]) && ("/" %in% strsplit(m[1,1],"")[[1]])) { #if alleles are merged to 1 genotype per column 
            splitvec = unlist(strsplit(m, "/", fixed=T))
            nrows = nrow(m)
            msplit = matrix(missing, ncol=2*n, nrow=nrows)
            msplit[, 2*seq_len(n)-1] = splitvec[2*seq_len(n*nrows) - 1]
            msplit[, 2*seq_len(n)] = splitvec[2*seq_len(n*nrows)]
            m = msplit
        }
        if(ncol(m)%%2 != 0) stop("Uneven number of marker allele columns")
        nMark = ncol(m)/2
        if(!is.null(annotations)) {
            if(length(annotations) == 2 && !is.null(names(annotations))) annotations = rep(list(annotations), nMark) # if given attrs for a single marker
            else if(length(annotations) != nMark) stop("Length of marker annotation list does not equal number of markers.")
            
            markerdata_list = lapply(1:nMark, function(i) {
                if(is.null(attribs <- annotations[[i]])) return(NULL)
                mi = m[, c(2*i-1, 2*i), drop=FALSE]
                do.call(.createMarkerObject, c(list(matr=mi, missing=missing), attribs))
            })
        }
        else {
            markerdata_list = lapply(1:nMark, function(i) {
                mi = m[, c(2*i-1, 2*i), drop=FALSE]
                .createMarkerObject(mi, missing=missing)
            })
        }
        markerdata_list[sapply(markerdata_list, is.null)] = NULL
        class(markerdata_list) = "markerdata"
    }
    x$nMark = length(markerdata_list) 
    x$markerdata = markerdata_list
    x$available = if (x$nMark>0) x$orig.ids[rowSums(do.call(cbind, markerdata_list)) > 0] else numeric(0)
    x
}

.createMarkerObject = function(matr, missing, alleles=NULL, afreq=NULL, chrom=NA, pos=NA, name=NA) {
    if(is.null(alleles)) {
        vec = as.vector(matr)
        alleles = unique.default(vec[vec != missing])
        if (length(alleles)==0) alleles = 1
    }
    if(!is.numeric(alleles) && !any(grepl("[^0-9]", alleles)))
        alleles = as.numeric(alleles)
    all_ord = order(alleles)
    alleles = alleles[all_ord]
    nalleles = length(alleles)
    if(is.null(afreq)) afreq = rep.int(1, nalleles)/nalleles 
    else {
        if(length(afreq) != nalleles) stop("Number of alleles don't match length of frequency vector")
        if(round(sum(afreq),2) != 1) warning(paste("Allele frequencies for marker", name," do not sum to 1:", paste(afreq, collapse=", ")))
        afreq = afreq[all_ord]
    }
    m_obj = match(matr, alleles, nomatch=0)
    attributes(m_obj) = list(dim=dim(matr), missing=missing, alleles=as.character(alleles), nalleles=nalleles, afreq = afreq, chrom=chrom, pos=pos, name=name, class="marker") 
    m_obj
}


.readMap = function(map, dat, freq, verbose, numerical=FALSE) { #TODO: numerical not used
    stopifnot(!is.null(map), !is.null(dat))
    if(is.character(map) && length(map)==1) {
        rawmap = read.table(map, header=FALSE, comment.char="", colClasses="character")
        if(!any(1:9 %in% strsplit(rawmap[[1]][1],"")[[1]])) rawmap = rawmap[-1,] #If no number occurs in first entry, first row is assumed to be header.
    }
    else  
        rawmap = map 
    rawmap[[1]][rawmap[[1]]=="X"] = 23
    rawmap[[1]][rawmap[[1]]=="Y"] = 24
    rawmap[[1]][!rawmap[[1]] %in% 1:24] = NA
    mapnames = as.character(rawmap[[2]])
    map1 = matrix(as.numeric(c(rawmap[[1]], rawmap[[3]])), ncol=2, dimnames = list(mapnames, c('CHR','POS'))) 

    if(is.character(dat) && length(dat)==1) #if a file name
        rawdat = read.table(dat, header=FALSE, comment.char="", colClasses="character")
    else  
        rawdat = dat
    datnames = as.character(rawdat[rawdat[[1]] == "M", 2])  #names of all markers in map
    
    if(is.character(freq) && length(freq)==1) {
        rawfreq = as.matrix(read.table(freq, header=FALSE, colClasses="character", fill=T))
        NROW = nrow(rawfreq)
        if(rawfreq[2,1]=="F") {
            freqnames = rawfreq[seq(1, NROW, by=2), 2]
            freqmatr = matrix(as.numeric(rawfreq[seq(2, NROW, by=2), -1]), nrow=NROW/2, dimnames=list(freqnames, NULL))
            isna = is.na(freqmatr)
            nalls = rowSums(!isna)
            allelmatr = col(freqmatr)
            allelmatr[isna] = NA
        }
        else { #if rawfreq[2,1]=="A" ... i.e. long format! See MERLIN tutorial about input files.
            Mrow = which(rawfreq[, 1]=="M")
            freqnames = rawfreq[Mrow, 2]

            nalls = c(Mrow[-1], NROW+1) - Mrow -1
            nM = length(Mrow)

            freqmatr = matrix(integer(), ncol=max(nalls), nrow=nM, dimnames=list(freqnames, NULL))
            indexmatr = cbind(rep(seq_len(nM), times=nalls), unlist(lapply(nalls, seq_len)))
            freqmatr[indexmatr] = as.numeric(rawfreq[-Mrow, 3])

            allalleles = rawfreq[-Mrow, 2, drop=F]
            numerical = !any(is.na(suppressWarnings(all_num <- as.numeric(allalleles))))
            if(numerical) {
                allelmatr = matrix(numeric(), ncol=max(nalls), nrow=nM, dimnames=list(freqnames, NULL))
                allelmatr[indexmatr] = all_num
            } else {
                allelmatr = matrix(character(), ncol=max(nalls), nrow=nM, dimnames=list(freqnames, NULL))
                allelmatr[indexmatr] = allalleles
            }
        }
        seqlist = lapply(1:max(nalls), seq_len) # precomputing vectors to speed up mapinfo further down
    }
    else if(is.list(freq))    {
        stop("Unexpected frequency object in readmap()")
        freqlist = freq #old
    }
    else freqmatr = NULL
    
    mapMatch = match(datnames, mapnames, nomatch=0)
    if(verbose && any(NAs <- mapMatch==0))
        cat("Deleting the following marker(s), which are not found in the map file:\n", paste(datnames[NAs], collapse="\n"), "\n", sep="")
    
    if(is.null(freqmatr))
        annotations = lapply(seq_along(datnames), function(i) {
            if((mapmatch = mapMatch[i]) == 0) return(NULL)
            list(chrom = map1[mapmatch,1], pos = map1[mapmatch,2], name=datnames[i])
        })   
    else {
        freqMatch = match(datnames, freqnames, nomatch=0)
        annotations = lapply(seq_along(datnames), function(i) {
            if((mapmatch <- mapMatch[i]) ==0) return(NULL)
            if((fmatch <- freqMatch[i]) ==0) alleles = afreq = NULL
            else {
                sq = seqlist[[nalls[fmatch]]]
                alleles = allelmatr[fmatch, sq] ###TODO: CHECK THIS!!
                afreq = as.vector(freqmatr[fmatch, sq])
            }
            list(alleles=alleles, afreq=afreq, chrom = map1[mapmatch,1], pos = map1[mapmatch,2], name=datnames[i])
        })
   }
   annotations
}


marker = function(x, ..., allelematrix, alleles=NULL, afreq=NULL, missing=0, chrom=NA, name=NA, pos=NA) {
    arglist = list(...)
    n = length(arglist)
    if(n == 0) {
        if(missing(allelematrix)) m = matrix(missing, ncol=2, nrow=x$nInd)
        else {
            m = as.matrix(allelematrix)
            stopifnot(nrow(m)==x$nInd, ncol(m)==2)
        }
    }
    if(n > 0) {
        if(!missing(allelematrix)) stop("Syntax error. See ?marker.")
        if(n == 1 && length(arglist[[1]]) > 2) stop("Syntax error. See ?marker.")
        if(n > 1 && n%%2 != 0) stop("Wrong number of arguments.")
        
        fill = {if(n==1) arglist[[1]] else missing}
        m = matrix(fill, ncol=2, nrow=x$nInd, byrow=TRUE)  #create marker matrix with all individuals equal.
        
        for(i in (seq_len(n/2)*2 - 1)) {   # odd numbers up to n-1.
            ids = arglist[[i]]
            geno = arglist[[i+1]]
            for(j in match(ids, x$orig.ids)) m[j, ] = geno
        }
    }
    
    .createMarkerObject(m, missing=missing, alleles=alleles, afreq=afreq, chrom=chrom, name=name, pos=pos)
}

removeMarkers <- function(x, markers=NULL, markernames=NULL, chroms=NULL, from=NULL, to=NULL) {
    if(is.null(markers)) markers = getMarkers(x, markernames, chroms, from, to)
    if(is.null(markers) || length(markers)==0) return(x)
    m = x$markerdata
    m[markers] = NULL
    setMarkers(x, m)
}

.prettyMarkers = function(m, alleles=NULL, sep="", missing=NULL, singleCol=FALSE, sex) {
    if(is.null(m)) return(m) 
    if(is.matrix(m)) m = list(m)
    if((n <- length(m))==0) return(m)
    if(is.null(alleles)) alleles = lapply(m, attr, 'alleles')
    else {
        if(!is.atomic(alleles)) stop("The parameter 'alleles' must be NULL, or an atomic vector.")
        if(length(alleles) < max(unlist(lapply(m, attr, 'nalleles')))) stop("The indicated 'alleles' vector has too few alleles for some markers.")
        alleles = rep(list(alleles), n)
    }
    if(is.null(missing)) missing = unlist(lapply(m, attr, 'missing')) 
    else {
        if(!is.atomic(missing) || length(missing) != 1) stop("The parameter 'mising' must be NULL, or a numeric/character of length 1.") 
        missing = rep(missing, n)
    }
    mNames = unlist(lapply(m, attr, 'name')); mNames[is.na(mNames)]=""
    pretty_m = do.call(c, lapply(seq_len(n), function(i)  c(missing[i], alleles[[i]])[m[[i]]+1]))
    dim(pretty_m) = c(length(pretty_m)/(2*n), 2*n)
    if (singleCol) {
        al1 = pretty_m[, 2*seq_len(n)-1, drop=FALSE]
        al2 = pretty_m[, 2*seq_len(n), drop=FALSE]
        m.matrix = matrix(paste(al1, al2, sep=sep), ncol=n )
        chrom_X = unlist(lapply(m, function(mm) identical(23L, as.integer(attr(mm, 'chrom')))))
        if(any(chrom_X)) {
            males = (sex==1)
            if (!all(hh <- al1[males, chrom_X] == al2[males, chrom_X])) warning("Male heterozygosity at X-linked marker detected.")
            m.matrix[males, chrom_X] = al1[males, chrom_X]
        }
        colnames(m.matrix) = mNames
        return(m.matrix)
    }
    else {
        nam = rep(mNames, each=2); nam[nzchar(nam)] = paste(nam[nzchar(nam)], 1:2, sep="_")
        colnames(pretty_m) = nam
        return(pretty_m)
    }
}

    
modifyMarker = function(x, marker, ids, genotype, alleles, afreq, chrom, name, pos) {
    if(inherits(marker, "marker")) {
        if(nrow(marker) != x$nInd) stop("Wrong dimensions of marker matrix.")
        m = marker
    }
    else {
        if (!is.numeric(marker) || length(marker) != 1) stop("Argument 'marker' must be a single integer or an object of class 'marker'.")
        if (marker > x$nMark) stop("Indicated marker does not exist.")
        m = x$markerdata[[marker]]
    }
    mis = attr(m, 'missing')
    
    if(!missing(alleles)) {
        stopifnot(is.atomic(alleles), is.numeric(alleles) || is.character(alleles))
        #if(is.numeric(alleles)) alleles = as.character(alleles)
        if(attr(m, 'missing') %in% alleles) stop("The 'missing allele' character cannot be one of the alleles.")
        lena = length(alleles)
        if(lena == attr(m, 'nalleles'))
            attr(m, 'alleles') = as.character(alleles)
        else {
            num_als = unique.default(as.vector(m[m!=0]))
            if(lena < length(num_als)) stop("Too few alleles.")
            pm = matrix(c(0, attr(m, 'alleles'))[m + 1], ncol=2)
            m = .createMarkerObject(pm, missing=0, alleles=alleles, afreq=rep(1, lena)/lena, chrom=attr(m, 'chrom'), pos=attr(m, 'pos'), name=attr(m, 'name'))
        }
    }
    if(!missing(afreq)) {
        if(round(sum(afreq),2)!=1) stop("The allele frequencies don't sum to 1.")
        if(length(afreq) != attr(m, 'nalleles')) stop("The length of allele frequency vector doesn't equal the number of alleles.")
        if(is.null(names(afreq))) 
            attr(m, 'afreq') = afreq
        else
            if (all(names(afreq) %in% attr(m, 'alleles'))) attr(m, 'afreq') = afreq[attr(m, 'alleles')] else stop("The names of the frequency vector don't match the allele names.")
    }
    
    changegeno = sum(!missing(ids), !missing(genotype))
    if(changegeno == 1) stop("The parameters 'ids' and 'genotype' must either both be NULL or both non-NULL.")
    if(changegeno == 2) {
        if(!is.atomic(genotype) || length(genotype)>2) stop("The 'genotype' parameter must be a numeric or character vector of length 1 or 2.")
        pm = .prettyMarkers(list(m))
        for (i in .internalID(x, ids)) pm[i, ] = genotype
        if(!all(as.vector(pm) %in% c(attr(m, 'alleles'), mis))) stop("Unknown allele(s). Please specify allele names using the 'alleles' argument.")
        attribs = attributes(m)
        m = match(pm, attr(m, 'alleles'), nomatch=0)
        attributes(m) = attribs
    }
    if(!missing(chrom)) attr(m, 'chrom') = chrom
    if(!missing(name)) attr(m, 'name') = name
    if(!missing(pos)) attr(m, 'pos') = pos
    if(inherits(marker, "marker")) return(m)
    else {
        x$markerdata[[marker]] = m
        x
    }
}

getMarkers = function(x, markernames=NULL, chroms=NULL, from=NULL, to=NULL) {
    mnos = seq_len(x$nMark)
    if(!is.null(markernames)) {
        if(length(markernames)==0) return(numeric(0))
        mnos = mnos[match(markernames, unlist(lapply(x$markerdata, function(m) attr(m, 'name'))), nomatch=0)]
    }
    if(!is.null(chroms)) {
        if(length(chroms)==0) return(numeric(0))
        mnos = mnos[unlist(lapply(x$markerdata[mnos], function(m) attr(m, 'chrom'))) %in% chroms]
    }
    if(!is.null(from)) 
        mnos = mnos[unlist(lapply(x$markerdata[mnos], function(m) attr(m, 'pos'))) >= from]
    if(!is.null(to)) 
        mnos = mnos[unlist(lapply(x$markerdata[mnos], function(m) attr(m, 'pos'))) <= to]
    mnos
}

swapGenotypes = function(x, ids){
    stopifnot(length(ids)==2)
    ids = .internalID(x, ids)
    y = as.matrix(x)
    y[ids, -(1:6)] = y[ids[2:1], -(1:6)]
    restore_linkdat(y)
}

modifyMarkerMatrix = function(x, ids, new.alleles){
    ids = .internalID(x, ids)
    y = as.matrix(x)
    y[ids, -(1:6)] = new.alleles
    restore_linkdat(y)
}


mendelianCheck = function(x, remove=FALSE, verbose=!remove) {
   if(inherits(x, 'singleton')) return(numeric(0))
    
    trioCheckFast = function(fa, mo, of) {
      even = 2*seq_len(length(fa)/2); odd = even - 1
        fa_odd = fa[odd]; fa_even = fa[even]; mo_odd = mo[odd]; mo_even = mo[even]
        of_odd = of[odd]; of_even = of[even];
        fa0 = (fa_odd == 0 | fa_even == 0)
        mo0 = (mo_odd == 0 | mo_even == 0)
        of_odd0 = (of_odd == 0); of_even0 = (of_even == 0)
        ff1 = (fa0 | of_odd0 | of_odd == fa_odd | of_odd == fa_even)
        ff2 = (fa0 | of_even0 | of_even == fa_odd | of_even == fa_even)
        mm1 = (mo0 | of_odd0 | of_odd == mo_odd | of_odd == mo_even)
        mm2 = (mo0 | of_even0 | of_even == mo_odd | of_even == mo_even)
        (ff1 & mm2) | (ff2 & mm1)
    }
    
    maleXHomoz = function(of)     {
        even = 2*seq_len(length(of)/2); odd = even - 1
        of[odd] == of[even]
    }
    
    maleXCheck = function(mo, of)    {
        even = 2*seq_len(length(of)/2); odd = even - 1
        mo_odd = mo[odd]; mo_even = mo[even]
        of = of[odd]
        mo0 = (mo_odd == 0 | mo_even == 0)
        of == 0 | mo0 | of == mo_odd | of == mo_even
    }
    
    nuclearAllelCheck = function(parents, offs) { #both matrices with 2*N columns
        even = 2*seq_len(ncol(parents)/2)
        unlist(lapply(even, function(i) {
            offs_als = unique.default(offs[,(i-1):i])
            offs_count = length(offs_als[offs_als > 0])
    
            par_als = parents[,(i-1):i]
            iszero = (par_als == 0); nzero = sum(iszero)
            if(nzero==0)     parents_count = length(unique.default(par_als))
            else { 
                par_als = unique.default(par_als[!iszero])
                parents_count = length(.myintersect(par_als, offs_als)) + nzero
            }
            offs_count <= parents_count
        }))
    }    
    
    ped = x$pedigree
    parents = unique(ped[, 2:3]) 
    parents = parents[ -match(0, parents[,1]), , drop=FALSE]
    subnucs = lapply(nrow(parents):1, function(i) { 
        par = parents[i,];  
        c(fa=par[[1]], mo=par[[2]], offs=as.vector(ped[,1])[ which(ped[,2]==par[[1]] & ped[,3]==par[[2]], useNames=FALSE) ]) 
    })
    
    chromX = unlist(lapply(x$markerdata, function(mm) identical(23L, as.integer(attr(mm, 'chrom')))))
    which_AUT = which(!chromX)
    which_X = which(chromX)
   
    errorlist = rep(list(numeric(0)), nrow(ped))
    names(errorlist) = x$orig.ids
    nuc_errors = numeric() # container for allele count errors...belongs to the whole subnuc.
    
    ### AUTOSOMAL
    if(length(which_AUT) > 0) {
        if (verbose) cat("\n### Checking autosomal markers ###\n")
        mdat = do.call(cbind, x$markerdata[!chromX])
        for (sub in subnucs) {
            fa = mdat[sub[[1]], ]; mo = mdat[sub[[2]], ]; offs = sub[-(1:2)]
            for (of in offs) {
                new_errors = which_AUT[!trioCheckFast(fa, mo, mdat[of, ])]
                if(length(new_errors) > 0) {
                    errorlist[[of]] = c(errorlist[[of]], new_errors)
                    if(verbose) cat("Individual ", x$orig.ids[of], " incompatible with parents for ", length(new_errors), " markers: ",
                               paste(new_errors, collapse=", "),"\n", sep="")
                }
            }
            if (length(offs) > 1) {
                new_errors = which_AUT[!nuclearAllelCheck(parents = mdat[sub[1:2], ], offs = mdat[sub[-(1:2)], ])]
                if(length(new_errors) > 0) {
                    nuc_errors = c(nuc_errors, new_errors)
                    if(verbose) cat("Offspring of ", x$orig.ids[sub[1]], " and ", x$orig.ids[sub[2]], " have too many alleles for ", length(new_errors), " markers: ", 
                            paste(new_errors, collapse=", "),"\n", sep="")
                }
            }
        }
    }
   
    ### X
    if(length(which_X) > 0) {
        if (verbose) cat("\n### Checking markers on the X chromosome ###\n")
        sex = x$pedigree[, 'SEX']
        mdat = do.call(cbind, x$markerdata[chromX])
        for (sub in subnucs) {
            fa = mdat[sub[[1]], ]; mo = mdat[sub[[2]], ]; offs = sub[-(1:2)]
            for (of in offs) {
                ofdat = mdat[of, ]
                if(sex[of]==1) {
                    new_errors = which_X[!maleXHomoz(ofdat)]
                    if(length(new_errors) > 0) {
                        errorlist[[of]] = c(errorlist[[of]], new_errors)
                        if(verbose) cat("Male ", x$orig.ids[of], " heterozygous for ", length(new_errors), " X-linked markers: ", paste(new_errors, collapse=", "), "\n", sep="")
                    }
                    new_errors = which_X[!maleXCheck(mo, ofdat)]
                    if(length(new_errors) > 0) {
                        errorlist[[of]] = c(errorlist[[of]], new_errors)
                        if(verbose) cat("Male ", x$orig.ids[of], " incompatible with mother for ", length(new_errors), " X-linked markers: ", paste(new_errors, collapse=", "), "\n", sep="")
                    }
                }
                else {
                    new_errors = which_X[!trioCheckFast(fa, mo, ofdat)]
                    if(length(new_errors) > 0) {
                        errorlist[[of]] = c(errorlist[[of]], new_errors)
                        if(verbose) cat("Female ", x$orig.ids[of], " incompatible with parents for ", length(new_errors), " X-linked markers: ", paste(new_errors, collapse=", "), "\n", sep="")
                    }
                }
            }
            if (length(offs) > 2) {
                new_errors = which_X[!nuclearAllelCheck(parents = mdat[sub[1:2], ], offs = mdat[offs, ])] #TODO: fix this for X-linked
                if(length(new_errors) > 0) {
                    nuc_errors = c(nuc_errors, new_errors)
                    if(verbose) cat("Offspring of ", x$orig.ids[sub[1]], " and ", x$orig.ids[sub[2]], " have too many alleles for ", length(new_errors), " markers: ", 
                            paste(new_errors, collapse=", "),"\n", sep="")
                }
            }
        }
    }
    err_index = sort.int(unique.default(c(unlist(errorlist), nuc_errors)))
   
    if(remove) return(removeMarkers(x, err_index))
    else return(err_index)
}

merlinUnlikely <- function(x, remove=FALSE, verbose=!remove) {
    merlin(x, model=F, options="--error --prefix _merlinerror", verbose=verbose)
    err = read.table("_merlinerror.err", header=T, as.is=T)
    unlink("_merlinerror.err")
    if(verbose) print(err)
    if(remove) return(removeMarkers(x, markernames=err$MARKER))
    else return(getMarkers(x, markernames=err$MARKER))
}

.merlin.unlikely = merlinUnlikely