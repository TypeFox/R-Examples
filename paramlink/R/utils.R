.mysetdiff = function(x,y) x[match(x,y,0L)==0L]
.myintersect = function(x,y) y[match(x, y, 0L)]

.internalID = function(x, orig.ids) {
    internal_ids = match(orig.ids, x$orig.ids)
    if(any(is.na(internal_ids))) stop(paste("Indicated ID(s) not among original ID labels:", paste(orig.ids[is.na(internal_ids)], collapse=",")))
    internal_ids
}

.getSex = function(x, orig.ids) as.vector(x$pedigree[.internalID(x, orig.ids), 'SEX'])

.comb2 = function(n) {
        if (n<2) return(matrix(nrow=0, ncol=2))
        v1 = rep.int(seq_len(n-1), (n-1):1)
        v2 = NULL
        for (i in 2:n) v2 = c(v2, i:n)
        cbind(v1,v2, deparse.level=0)
} 

allGenotypes = function(n)
    rbind(cbind(seq_len(n), seq_len(n)), .comb2(n))

.rand01 = function(n) sample.int(2, size=n, replace=T) - 1 #random 0/1 vector of length n.
    
.prettycat = function(v, andor)
    switch(min(len <- length(v), 3), toString(v), paste(v, collapse=" and "), paste(paste(v[-len], collapse=", "), andor, v[len]))
    
fast.grid = function(argslist, as.list=FALSE) {
    nargs <- length(argslist) 
    orep <- nr <- prod(lengths(argslist))
    if (nargs == 0L || nr == 0L) return(matrix(ncol=0, nrow=0))
    
    rep.fac <- 1L
    res <- NULL
    for (x in argslist) {
        nx <- length(x)
        orep <- orep/nx
        res <- c(res, x[rep.int(rep.int(seq_len(nx), rep.int(rep.fac, nx)), orep)]) #this is res[, i]
        rep.fac <- rep.fac * nx
    }
    dim(res) <- c(nr, nargs)
    if(as.list) res = lapply(seq_len(nr), function(r) res[r,])
    res
}
.my.grid = fast.grid

offspring = function(x, id, original.id=TRUE) {
    if(original.id) id = .internalID(x, id)
    p = x$pedigree
    offs_rows = p[, 1 + p[id, 'SEX'] ] == id
    if(original.id) x$orig.ids[offs_rows] else (1:x$nInd)[offs_rows]
}    


spouses = function(x, id, original.id=TRUE) { #returns a vector containing all individuals sharing offspring with 'i'.
    internal_id = ifelse(original.id, .internalID(x, id), id)
    p = x$pedigree
    offs_rows = p[, 1 + p[internal_id, 'SEX'] ] == internal_id
    spou = unique.default(p[offs_rows, 4 - p[internal_id, 'SEX']])  #sex=1 -> column 3; sex=2 -> column 2.
    if (original.id) return(x$orig.ids[spou]) else return(spou)
}

related.pairs = function(x, relation=c('parents', 'siblings', 'grandparents', 'nephews_nieces', 'cousins', 'spouses', 'unrelated'), 
                         available=F, original.id=T, interfam=c("none", "founders", "all"), ...) {
    relation = match.arg(relation)
    interfam = match.arg(interfam)
    func = function(...) get(relation)(...)
    
    if(!inherits(x, 'linkdat') && length(x)>0 && all(unlist(lapply(x, inherits, 'linkdat')))) {
        res = do.call(rbind, lapply(x, function(xx) {
            pairs = related.pairs(xx, relation, available, original.id, ...)
            if(nrow(pairs) > 0) {
                newpairs = paste(xx$famid, pairs, sep="-")
                dim(newpairs) = dim(pairs)
                newpairs
            }
            else pairs
        }))
        if(relation=="unrelated" && interfam != "none") {
            avail = lapply(x, function(xx) {
                ids = if(available) xx$available else xx$orig.ids
                if(interfam == "founders") ids = intersect(ids, xx$orig.ids[xx$founders])
                if(length(ids)==0) NULL else paste(xx$famid, ids, sep="-")})
            avail = avail[!sapply(avail, is.null)]
            fampairs = data.frame(t(.comb2(length(avail)))) # enable use of lapply below
            interfam = do.call(rbind, lapply(fampairs, function(p) fast.grid(avail[p])))
            res = rbind(res, interfam)
        }
        return(res)
    }
    
    res = NULL
    for(i in 1:x$nInd) {
        rels = func(x, i, original.id=F, ...)
        rels = rels[rels != i]
        res = rbind(res, cbind(rep.int(i, length(rels)), rels, deparse.level=0))
    }
    res[res[,1] > res[,2],] = res[res[,1] > res[,2], 2:1]
    res = unique(res)
    if(available) {
        avail = .internalID(x, x$available)
        res = res[res[,1] %in% avail & res[,2] %in% avail, , drop=F]
    }
    if(original.id) {
        res = matrix(x$orig.ids[res], ncol=2)
        res[res[,1] > res[,2],] = res[res[,1] > res[,2], 2:1]  # sort again
    }
    res[order(res[,1], res[,2]),]
}

unrelated = function(x, id, original.id=TRUE) {
    if(!original.id) id = x$orig.ids[id]
    ancs = c(id, ancestors(x, id))
    rel = unique.default(unlist(lapply(ancs, function(a) c(a, descendants(x, a, original.id=TRUE)))))
    unrel = setdiff(x$orig.ids, rel)
    if(!original.id) unrel = .internalID(x, unrel)
    unrel
}
    
    
parents = function(x, id, original.id = TRUE) {
    grandparents(x, id, degree=1, original.id=original.id)
}

grandparents = function(x, id, degree=2, original.id = TRUE) {
    if (original.id) id = .internalID(x, id)
    p = x$pedigree
    gp = id
    for(i in seq_len(degree))
        gp = p[gp, 2:3]
    if (original.id) x$orig.ids[gp]
    else (1:x$nInd)[gp]
}

siblings = function(x, id, half = NA, original.id = TRUE) {
    if (original.id) id = .internalID(x, id)
    p = x$pedigree
    fa = p[id, 'FID']; mo = p[id, 'MID']
    if(fa==0 && mo==0) return(numeric()) 
    samefather = p[, 'FID'] == fa
    samemother = p[, 'MID'] == mo
    sib_rows = if(is.na(half))   samefather | samemother
               else if(half) xor(samefather,  samemother)
               else              samefather & samemother
    sib_rows[id] = FALSE
    if (original.id) x$orig.ids[sib_rows]
    else (1:x$nInd)[sib_rows]
}

cousins = function(x, id, degree=1, removal=0, half = NA, original.id = TRUE) {
    if (original.id) id = .internalID(x, id)
    gp = grandparents(x, id, degree=degree, original.id=FALSE)
    uncles = unique.default(unlist(lapply(gp, siblings, x=x, half=half, original.id=FALSE)))
    cous = uncles
    for(i in seq_len(degree+removal)) 
        cous = unique.default(unlist(lapply(cous, offspring, x=x, original.id=FALSE)))
    if (original.id) cous = x$orig.ids[cous]
    cous
}

nephews_nieces = function(x, id, removal=1, half = NA, original.id = TRUE) {
    cousins(x, id, degree=0, removal=removal, half=half, original.id=original.id)
}

ancestors = function(x, id) { #climbs up the pedigree storing parents iteratively. (Not documented: Accepts id of length > 1)
    if (inherits(x, 'linkdat')) {
        p = x$pedigree
        orig_ids = x$orig.ids
        ids_int = .internalID(x, id)
    }
    else if(is.matrix(x) && c('ID', 'FID', 'MID') %in% colnames(x)) {
        p = x
        orig_ids = p[, 'ID']
        ids_int = match(id, orig_ids)
    }
    else stop("x must be either a linkdat object or a matrix whose colnames include 'ID', 'FID' and 'MID'")
    p = relabel(p, 1:nrow(p))
   
    ancest = numeric(0)
    up1 = as.numeric(p[ids_int, c('FID', 'MID')])
    up1 = up1[up1 > 0 & up1 <= nrow(p)]  #NB: Avoids pedigree errors without warning! Should be caught in .checkped anyway  
    up1 = up1[!duplicated.default(up1)]
    while(length(up1) > 0) {
        ancest = c(ancest, up1)
        up1 = .mysetdiff(as.numeric(p[up1, c('FID', 'MID')]), ancest)
    }
    ancest = sort.int(ancest[(ancest != 0) & !duplicated(ancest)])
    return(orig_ids[ancest])
}

descendants = function(x, id, original.id=TRUE) { 
    internal_id = ifelse(original.id, .internalID(x, id), id)
    nextgen <- desc <- offspring(x, internal_id, original.id=FALSE)
    while(TRUE) {
        nextgen <- unlist( lapply(nextgen, offspring, x=x, original.id=FALSE))
        if (length(nextgen)==0) break
        desc <- c(desc, nextgen)
    }
    desc = unique.default(sort.default(desc))
    if (original.id) return(x$orig.ids[desc]) else return(desc)
}

.generations = function(x) {#linkdat object
    max(vapply(unlist(.descentPaths(x, x$founders, original.ids=FALSE), recursive=F), length, 1))
}

.descentPaths = function(x, ids, original.ids = TRUE)  {
    if (original.ids) ids = .internalID(x, ids)
    offs = lapply(1:x$nInd, offspring, x=x, original.id=FALSE)
    lapply(ids, function(id) {
        res=list(id)
        while(TRUE) {
            newoffs = offs[sapply(res, function(path) path[length(path)])]
            if (length(unlist(newoffs))==0) break
            nextstep = lapply(1:length(res), function(r) if (length(newoffs[[r]])==0) res[r] else lapply(newoffs[[r]], function(kid) c(res[[r]], kid) ))
            res = unlist(nextstep, recursive=FALSE)
        }
        if (original.ids) lapply(res, function(internal_vec) x$orig.ids[internal_vec]) else res
    })
}

pedigreeLoops = function(x) {
    dls = .descentPaths(x, 1:x$nInd, original.ids=FALSE)
    loops = list()
    for (id in 1:x$nInd) {
        if (length(dl <- dls[[id]])==1) next
        pairs = .comb2(length(dl))
        for (p in 1:nrow(pairs)) {
            p1 = dl[[pairs[p,1]]]; p2 = dl[[pairs[p,2]]]
            if (p1[2]==p2[2]) next
            inters = p1[match(p2,p1,0L)][-1]  #intersecting path members, excluding id
            if (length(inters)==0) next
            else {
                top=x$orig.ids[p1[1]]; bottom=x$orig.ids[inters[1]]
                pathA = p1[ seq_len(which(p1==inters[1]) - 2) + 1 ] #without top and bottom. Seq_len to deal with the 1:0 problem.
                pathB = p2[ seq_len(which(p2==inters[1]) - 2) + 1 ]
                loops = c(loops, list(list(top=top, bottom=bottom, pathA = x$orig.ids[pathA], pathB=x$orig.ids[pathB])))
            }
        }
    }
    unique(loops)
}


breakLoops = function(x, loop_breakers) {
    if(inherits(x,'singleton')) stop("This function does not apply to singleton objects.")
   stopifnot(is.numeric(loop_breakers), length(loop_breakers)>0)
    if(any(loop_breakers %in% x$orig.ids[x$founders])) stop("Pedigree founders cannot be loop breakers.")

    pedm = as.matrix(x)#data.frame(x, famid=T, missing=0)
    attrs = attributes(pedm) #all attributes except 'dim'
    dup_pairs = x$loop_breakers #normally = NULL at this stage 
    for(id in loop_breakers) {
        dup_id = max(pedm[, 'ID']) + 1
        dup_pairs = rbind(dup_pairs, c(id, dup_id))
        intern = match(id, pedm[,'ID'])  #don't use .internalID here, since pedm changes all the time
        sex_col = pedm[intern, 'SEX'] + 2 # FID column if 'intern' is male; MID if female
        
        pedm = pedm[c(1:intern, intern, (intern+1):nrow(pedm)), ]
        pedm[intern+1, c('ID', 'FID','MID')] = c(dup_id, 0, 0)
        pedm[pedm[, sex_col] == id, sex_col] = dup_id
    }
    newx = restore_linkdat(pedm, attrs=attrs) 
    newx$loop_breakers = dup_pairs
    newx
}

tieLoops = function(x) {
    dups = x$loop_breakers
    if(is.null(dups) || nrow(dups)==0) {cat("No loops to tie\n"); return(x)}
    if(!all(dups %in% x$orig.ids)) stop("Something's wrong: Duplicated individuals no longer in pedigree.")
    pedm = as.matrix(x)
    attrs = attributes(pedm)

    origs = dups[,1]; copies = dups[,2]
    pedm = pedm[-match(copies, pedm[,'ID']), ]
    for(i in 1:length(origs)) {
        orig = origs[i]; copy = copies[i]
        sex = pedm[pedm[,'ID']==orig, 'SEX']
        pedm[pedm[, sex+2] == copy, sex+2] = orig 
    }
    restore_linkdat(pedm, attrs=attrs)
}

all.equal.linkdat = function(target, current, ...) {
    res = TRUE
    if(!all.equal(class(target), class(current))) {
        cat("Class attributes are not equal\n")
        res = FALSE
    }
    if(target$nMark != current$nMark) {
        cat("Unequal numbers of markers:", target$nMark, "vs.", current$nMark, "\n")
        res = FALSE
    }
    names.t = names(target)
    names.c = names(current)
    if(!setequal(names.t, names.c)) {
        if(length(first_miss <- setdiff(names.c, names.t)) > 0)
            cat("Missing slots in first object:", paste(first_miss, sep=", "), "\n")
        if(length(sec_miss <- setdiff(names.t, names.c)) > 0)
            cat("Missing slots in second object:", paste(sec_miss, sep=", "), "\n")
        res = FALSE
    }
    if(!setequal(target$orig.ids, current$orig.ids)) {
        cat("ID labels are not equal\n")
        res = FALSE
    }
   
    new_order = match(current$orig.ids, target$orig.ids)
    ped_targ = relabel(target$pedigree, target$orig.ids)[new_order, , drop=F]
    ped_curr = relabel(current$pedigree, current$orig.ids)
    if(!identical(ped_curr, ped_targ)) {
        cat("Pedigree topologies are not equal\n")
        res = FALSE
    }
   
    if(!setequal(target$available, current$available)) {
        cat("Unequal vectors of availability\n")
        res = FALSE
    }
    if(!res) return(res)
   
    if(target$nMark > 0) {
        mark_targ <- do.call(cbind, as.list(target$markerdata))[new_order, , drop=F]
        mark_curr <- do.call(cbind, as.list(current$markerdata))
        if(!isTRUE(all.equal(mark_targ, mark_curr))) {
            diffs = which(mark_targ != mark_curr, arr.ind=T)
            cat("Differences in the following markers:", sort(unique((diffs[,2]+1) %/% 2)), "\n")
            res = FALSE
        }
        markerattr_targ <- lapply(target$markerdata, attributes)
        markerattr_curr <- lapply(current$markerdata, attributes)
        if(!identical(markerattr_targ, markerattr_curr)) {
            diffattr = which(sapply(seq_along(markerattr_targ), function(i) !identical(markerattr_targ[[i]], markerattr_curr[[i]])))
            cat("Difference in marker attributes for marker", diffattr, "\n")
            res = FALSE
        }
    }
    return(res)
}


inbreeding = function(x) {
    ped = x$pedigree
    kin.coeff = kinship(id=ped[,'ID'], dadid=ped[,'FID'], momid=ped[,'MID'])
    inb.coeff = numeric()
    inb.coeff[x$founders] = 0
    inb.coeff[x$nonfounders] = sapply(x$nonfounders, function(i) kin.coeff[ped[i, 'FID'], ped[i, 'MID']])
    inb.coeff
}


#pedDistMatrix(x) computes a symmetrix integer matrix (d_ij), where d_ij is the pedigree distance (i.e., the length of the shortest pedigree path) between i and j. For example, the following relations have the indicated pedigree distances: parent/offspring    1, siblings    2, parents    2, uncle/niece    3, first cousins 4.
# .pedDistMatrix <- function(x) { 
    # d = matrix(0, ncol=x$nInd, nrow=x$nInd)
    # for (f in 1:x$nInd) {
        # d[ f, offspring(x, f, original.id=FALSE) ] <- 1
        # d[ f, x$pedigree[f, c('FID', 'MID')] ] <- 1
    # }
    # k=1
    # while (prod(d[col(d)!=row(d)])==0) { 
        # k=k+1
        # for (r in 1:nrow(cc <- .comb2(x$nInd))) {
            # a=cc[r,1]; b=cc[r,2]
            # if (d[a,b] > 0) next
            # if (any( d[a,]>0 & d[,b]>0 & d[a,]+d[,b] == k))
            # d[a,b] <- d[b,a] <- k
        # }
    # }
    # d
# }

# .getpaths = function(obj, i, j) {
    # if (inherits(obj, "linkdat")) d = .pedDistMatrix(obj)
    # else if (is.matrix(obj)) d = obj
    # else stop("'obj' must be either a matrix or a 'linkdat' object.")
    # if (i==j) return(i)
    # D = d[i,j]
    # res = matrix(i)
    # for (k in 1:D) {
        # tmp = list()
        # for (r in 1:nrow(res)) {
            # laststep = res[r, k]
            # newsteps = which( d[laststep, ] == 1 & d[, j] == D-k)
            # tmp[[r]] = cbind( matrix(res[r, 1:k], ncol=k, nrow=length(newsteps),byrow=TRUE), newsteps, deparse.level=0)
        # }
        # res = do.call(rbind, tmp)
    # }
    # res
# }