relabel <- function(x, new, old) {
    islinkdat = inherits(x,"linkdat")
    if(islinkdat) {
        if(length(new)==x$nInd && all(new==x$orig.ids)) return(x)
        ped = as.matrix(x)
        avail = attr(ped, "available")
    } 
    else ped = x
    
    orig.ids = ped[, 'ID']
    if(missing(old)) old = orig.ids
    stopifnot(is.numeric(old), is.numeric(new), length(old)==length(new), !0 %in% new, all(old %in% ped[, 'ID']))
    ped[match(old, orig.ids), 'ID'] = new
    
    parents = ped[, c('FID','MID')]
    ped[, c('FID','MID')][parents %in% old]   <- new[match(parents, old, nomatch=0)] #relabeling parents
    
    if(islinkdat) {
        oldavail = avail[avail %in% old]
        avail[avail %in% old] = new[match(oldavail, old)]
        attr(ped, "available") = avail
        return(restore_linkdat(ped))
    }
    else return(ped)
}

setAvailable = function(x, available) {
    x$available = sort(as.numeric(available))
    x
}

swapAvailable = function(x, ids) {
    ava = x$available
    new_ava = c(ava[!ava %in% ids], ids[!ids %in% ava])
    setAvailable(x, new_ava)
}
    
swapSex <- function(x, ids, verbose=TRUE) {
    ids = .internalID(x, ids)
    ids.spouses = unique(unlist(lapply(ids, spouses, x=x, original.id = FALSE)))
    if (!all(ids.spouses %in% ids)) {
        if(verbose) cat("Changing sex of the following spouses as well:", paste(x$orig.ids[setdiff(ids.spouses, ids)], collapse=", "), "\n")
        return(swapSex(x, x$orig.ids[union(ids, ids.spouses)]))
    }    
    pedm = as.matrix(x)
    pedm[ids, 'SEX'] <- (3 - pedm[ids, 'SEX'])
    offs = x$pedigree[,"FID"] %in% ids
    pedm[offs, c('FID', 'MID')] <- pedm[offs, c('MID', 'FID')]
    
    restore_linkdat(pedm)
}

swapAff <- function(x, ids, newval=NULL) {
    pedm = as.matrix(x)
    ids = .internalID(x, ids)
    if (is.null(newval))  newval <- (3 - pedm[ids, 'AFF'])
    
    pedm[ids, 'AFF'] <- newval
    restore_linkdat(pedm)
}


addOffspring <- function(x, father, mother, noffs, ids=NULL, sex=1, aff=1, verbose=TRUE) {
    p = as.matrix(x); attrs = attributes(p); nm = x$nMark
    taken <- oldids <- p[,'ID']
    if(!missing(father)) taken = c(taken, father)
    if(!missing(mother)) taken = c(taken, mother)
    if(!is.null(ids)) taken = c(taken, ids)
    max_id = max(taken)
    
    if(missing(father) && missing(mother)) stop("At least one parent must be an existing pedigree member.")
    if(missing(father))     father <- max_id <- max_id + 1 
    if(missing(mother))     mother <- max_id <- max_id + 1 
    if(any(!is.numeric(father), length(father)!=1))     stop("Argument 'father' must be a single integer.")
    if(any(!is.numeric(mother), length(mother)!=1))     stop("Argument 'mother' must be a single integer.")
    if(!any(c(father,mother) %in% oldids))    stop("At least one parent must be an existing pedigree member.")

    if (missing(noffs) && is.null(ids)) stop("Number of offspring not indicated.")
    if (missing(noffs)) noffs = length(ids)
    if (is.null(ids)) ids = (max_id+1):(max_id+noffs)
    if (length(ids)!=noffs) stop("Length of 'id' vector must equal number of offspring.")
    if (any(ids %in% oldids))    stop(paste("Individual(s)", ids[ids %in% oldids], "already exist(s)."))

    if(!father %in% oldids) {
        if(verbose) cat("Father: Creating new individual with ID", father, "\n")
        p = rbind(p, c(x$famid, father, 0, 0, 1, 1, rep.int(0, nm*2)))
    }
    if(!mother %in% oldids) {
        if(verbose) cat("Mother: Creating new individual with ID", mother, "\n")
        p = rbind(p, c(x$famid, mother, 0, 0, 2, 1, rep.int(0, nm*2)))
    }
    p = rbind(p, cbind(x$famid, ids, father, mother, sex, aff, matrix(0, ncol=nm*2, nrow=length(ids))))

    restore_linkdat(p, attrs=attrs)
}

addSon = function(x, parent, id=NULL, aff=1, verbose=TRUE) {
    if(.getSex(x, parent) == 1) 
        addOffspring(x, father=parent, noffs=1, sex=1, aff=aff, ids=id, verbose=verbose)
    else 
        addOffspring(x, mother=parent, noffs=1, sex=1, aff=aff, ids=id, verbose=verbose)    
}

addDaughter = function(x, parent, id=NULL, aff=1, verbose=TRUE) {
    if(.getSex(x, parent) == 1) 
        addOffspring(x, father=parent, noffs=1, sex=2, aff=aff, ids=id, verbose=verbose)
    else 
        addOffspring(x, mother=parent, noffs=1, sex=2, aff=aff, ids=id, verbose=verbose)    
}

addParents <- function(x, id, father, mother, verbose=TRUE) {
    if(length(id)>1) stop("Only one individual at the time, please")
    if(id %in% x$orig.ids[x$nonfounders]) stop(paste("Individual", id, "already has parents in the pedigree")) 
    
    p = as.matrix(x); attrs = attributes(p); nm = x$nMark
    oldids = p[,'ID']
    n = max(oldids)
    if (missing(father)) father <- n <- n+1    
    if (missing(mother)) mother = n + 1
    new.father = !father %in% oldids
    new.mother = !mother %in% oldids
    if(new.father && verbose) cat("Father: Creating new individual with ID", father, "\n")
    if(new.mother && verbose) cat("Mother: Creating new individual with ID", mother, "\n")
    
    int.id <- .internalID(x, id)
    p[int.id, c('FID', 'MID')] <- c(father, mother)

    if(new.father) 
        p = rbind(p, c(x$famid, father, 0, 0, 1, 1, rep.int(0, nm*2)))[append(1:nrow(p), nrow(p)+1, after=int.id - 1), ] #insert father before 'id'

    if(new.mother) 
        p = rbind(p, c(x$famid, mother, 0, 0, 2, 1, rep.int(0, nm*2)))[append(1:nrow(p), nrow(p)+1, after=int.id - 1 + as.numeric(new.father)), ] #insert mother before 'id'

    restore_linkdat(p, attrs=attrs)
}


removeIndividuals <- function(x, ids, verbose=TRUE) { #removes (one by one) individuals 'ids' and all their descendants. Spouse-founders are removed as well.
    if(any(!ids %in% x$orig.ids)) stop(paste("Non-existing individuals:", .prettycat(ids[!ids %in% x$orig.ids], "and")))  
    pedm = as.matrix(x)
    
    #founders without children after 'id' and 'desc' indivs are removed. The redundancy here does not matter.
    desc = numeric(0)
    for(id in ids) {
        desc = c(desc, dd <- descendants(x, id))
        if(verbose) cat("Removing", id, if(length(dd)>0) paste("and descendant(s):", .prettycat(dd, "and")), "\n")
    }

    leftover.spouses = setdiff(x$orig.ids[x$founders], c(ids, as.numeric(pedm[!x$orig.ids %in% c(ids, desc), c('FID','MID')])))  #founders that are not parents of remaining indivs
    if(verbose && length(leftover.spouses)>0) cat("Removing leftover spouse(s):", .prettycat(leftover.spouses, "and"), "\n")    

    remov = unique(c(ids, desc, leftover.spouses))
    restore_linkdat(pedm[-.internalID(x, remov), , drop=F], attrs=attributes(pedm))
}

branch <- function(x, id) {
   desc = descendants(x, id)
   spous = unlist(lapply(c(id, desc),spouses,x=x))
   subset(x, subset=c(id, desc, spous))
}

trim <- function(x, keep=c("available", "affected"), return.ids=FALSE, verbose=TRUE){
    keep = match.arg(keep)
    if(verbose) 
        cat("Trimming pedigree, keeping", keep, "individuals.")
    if(inherits(x, 'singleton') | (keep == "available" && length(x$available)==length(x$orig.ids))) {
        if(verbose) cat(" Removed: None\n")
        return(x)
    }
    mysetdiff = function(x,y) x[match(x,y,0L)==0L]
    
    y = linkdat(relabel(x$pedigree, x$orig.ids), verbose=F) # make a copy of x$ped, with original IDs
    y$available = x$available
    while(TRUE) {
        p = y$pedigree
        leaves = mysetdiff(p[, 'ID'], p[, c('FID','MID')])
        throw = switch(keep,
            available = mysetdiff(leaves, .internalID(y, y$available)), #nonavailable leaves
            affected = leaves[p[leaves, 'AFF'] != 2] #nonaffected leaves
        )
        if (length(throw)==0) break
        y = removeIndividuals(y, y$orig.ids[throw], verbose=FALSE)
    }

    remov = setdiff(x$orig.ids, y$orig.ids)
    if(return.ids) return(remov)
        
    store = as.matrix(x)
    trimmed = store[!store[, 'ID'] %in% remov, ]
    if(verbose) 
        cat(" Removed:", if(length(remov)>0) .prettycat(remov, "and") else "None", "\n")

    restore_linkdat(trimmed, attrs=attributes(store))
}

.merge.linkdat = function(x) {# list of linkdats
    if(!is.list(x) || !all(vapply(x, function(xx) inherits(xx, 'linkdat'), logical(1)))) stop("Input must be a list of linkdat objects")
    if(length(x)==1) return(x)
    mnames = lapply(x, function(xx) unlist(lapply(xx$markerdata, attr, 'name')))
    common = Reduce(intersect, mnames)
    lapply(x, function(xx) setMarkers(xx, xx$markerdata[getMarkers(xx, common)])) 
}

