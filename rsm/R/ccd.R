### Code for CCDs and related designs
# Modified by RVL starting Feb 10, 2011

### Find list of row indexes for each block in a design
.block.indices = function(design) {
    rsd = attr(design, "rsdes")
    if (!is.null(blknm <- rsd$block)) {
        blk = design[[blknm[1]]]
        if (length(blknm) > 1) for (nm in blknm[-1])
            blk = paste(blk, design[[nm]], sep=".")
        # Even now, it's possible to have a null blk, if that var not in data yet
        if (!is.null(blk)) 
            i.init = split(1:nrow(design), blk)
        else
            i.init = list(1:nrow(design))
    }
    else
        i.init = list(1:nrow(design))
    i.init
}

### Randomizes; adds columns "std.order" and "run.order" if it's a not-previously-randomized design
.randomize = function(design, randomize=TRUE) {
    OneToN = 1:nrow(design)
    if (is.na(match("std.order", names(design)))) {
        design$run.order = design$std.order = OneToN
        extcols = match(c("run.order","std.order"), names(design))
        oldcols = (1:ncol(design))[-extcols]
        design = `[.coded.data`(design, , c(extcols, oldcols))
    }
    if (randomize) {
        i.init = .block.indices(design)
        for (idx in i.init) {
            design[sample(idx), ] = design[idx, ]
            design$run.order[idx] = 1:length(idx)
        }
        row.names(design) = OneToN
    }
    design
}

### Calculate design moments needed for determining orthogonal blocks
### Computes (N, [i], [ii], [ij]) for each block (really just the [i1])
### Returns the matrix of [ii] values (col names are "N.obs" and pri vars, rows are blocks)
### For any block where any [i] != 0 or [ij] != 0, we return NA
### For orth blocks, need [ii]/N.obs constant on each block for each i
### (but not necessary for [ii] == [jj])
.orth.moments = function(design) {
    blkidx = .block.indices(design)
    if (is.null(rsd <- attr(design, "rsdes")))
        pricols = sapply(design, is.numeric)
    else
        pricols = match(rsd$primary, names(design))
    t(sapply(blkidx, function(idx) {
        N.obs = length(idx)
        X = as.matrix(design[idx, pricols])
        x1 = X[,1]
        moms = apply(X, 2, function(x) 
            c(i = sum(x), ii = sum(x^2), i1 = sum(x * x1)))
        if (any(abs(moms[1, ]) > 0.01) || (any(abs(moms[3,-1]) > .01)))
            moms[2, ] = NA
        c(N.obs=N.obs, moms[2, ])
    }))
}
    

### Calculate design moments needed for determining rotatability
### Computes design moments of order 4 or less and checks whether rotatability 
### is possible. We need all except [iiii] and [iijj] to be zero, 
### and the latter to be the same for all i <> j
### Returns iiii and iijj
### NOTE: For the design to be rotatable, we need [iiii]/[iijj] = 3
### For mixed moments, we actually check only [1j], [12k], etc.
.rot.moments = function(design) {
    if (is.null(rsd <- attr(design, "rsdes")))
        pricols = sapply(design, is.numeric)
    else
        pricols = match(rsd$primary, names(design))
    k = length(pricols)
    if (k < 2)
        stop("Need at least 2 primary variables for a rotatable design")
    X = as.matrix(design[, pricols])
    x1 = X[,1]; 
    x2 = X[,2]; 
    x3 = if (k==2) 0 else X[,3]
    moms = apply(X, 2, function(x) 
        c(i = sum(x), ij = sum(x*x1), 
          iii = sum(x^3), iij = sum(x^2*x1), ijk = sum(x*x1*x2),
          iiij = sum(x^3*x1), iijk = sum(x^2*x1*x2), ijkl = sum(x*x1*x2*x3),
          iiii = sum(x^4), iijj = sum(x^2 * x1^2)))
    # We need to check abs values of 1st 8 rows of moms
    moms[1:8,] = abs(moms[1:8,])
    # Zero out mixed parts that really aren't mixed
    moms["ij",1] = moms["iij",1] = moms["ijk",1] = moms["ijk",2] =
        moms["iiij",1] = moms["iijk",1] = moms["iijk",2] = 0
    if (k > 2) moms["ijkl",1:3] = 0
    # for checking purposes...
    moms["iijj", 1] = moms["iijj", 2]
    nzflag = any(moms[1:8,] > 0.01)
    if (nzflag || (sd(moms["iiii",]) > .1) || (sd(moms["iijj",]) > .1))
        moms["iiii", ] = moms["iijj", ] = NA
    c(iiii=moms["iiii",1], iijj=moms["iijj",2])
}

### Cube (1st-order) design -- this is a 2^k or 2^{k-p}, possibly replicated, plus center points
# inscribed may be boolean or anticipated value of alpha. If TRUE, we use sqrt(k) (spherical)
cube = function(basis, generators, n0=4, reps=1, coding, randomize=TRUE, 
    blockgen, bid=1, inscribed=FALSE)
{
    if (inherits(basis, "formula"))
        xvars = all.vars(basis[[length(basis)]])
    else if (is.numeric(basis))
        xvars = paste("x", 1:basis, sep="")
    else
        stop("'basis' must be an integer or a formula")
    
    CALL = match.call()
        
    args = lapply(xvars, function(nm) c(-1,1))
    names(args) = xvars
    des = do.call(expand.grid, args)
    
    if (!missing(generators)) {
        if (!is.list(generators)) generators = list(generators)
        for (gen in generators) {
            gen = as.character(gen)        
            des[[gen[[2]]]] = with(des, eval (parse (text = as.character(gen[[3]]))))
        }
    }
    
    k = ncol(des)
    xvars = names (des)
    # within-block reps
    if (reps > 1) {
        des = des[rep(1:nrow(des), each=reps), ]
    }
    
    rsd = list(primary=xvars)
    
    if (!missing(blockgen)) {
        blk = NULL
        if (!is.vector(blockgen)) {
            if (inherits(blockgen, "formula")) {
                bg = with(des, eval(blockgen[[length(blockgen)]])) ###blockgen = c(blockgen)
                if (is.numeric(bg)) blk = matrix(bg, nrow(des))                
            }
            else {
                bg = with(des, eval(blockgen))
                # should either give a list of formulas, or results
                if (is.numeric(bg)) 
                    blk = matrix(bg, nrow(des))
                else
                    blockgen = bg
            }
        }
        if (is.vector(blockgen) || is.null(blk))
            blk = with(des, sapply(blockgen, function(f) {
                if (is.character(f)) 
                    f = as.formula(paste("~", f))
                if (!inherits(f, "formula")) 
                    stop("block generators must be formulas or strings")
                if (length(f)==3) f = f[-2]
                eval(f[[2]])
            }))
        if (is.matrix(blk)) {
            fac = 2^(1:ncol(blk)) / 2   # = 1, 2, 4, ...
            blk = (blk %*% fac + 2^(ncol(blk)) + 1) / 2
        }
        else 
            blk = (blk + 3) / 2
        
        rep.rows = sapply(1:max(blk), function(b) which(blk==b)[1])
        rsd$blk.signs = des[rep.rows, ]
        rsd$bid = bid
        des = des[blk == bid, ]
    }
        
    # center points
    if (n0 > 0) {
        zero = as.data.frame(matrix(0, nrow=n0, ncol=k))
        names(zero) = names(des)
        des = rbind(des, zero)
    }
    
    if (inscribed) {
        if (is.logical(inscribed)) inscribed = sqrt(k)
        rsd$divisor = inscribed        
# Note: divisor is used only by star()
# Some trickery here. ccd() will give a negative value; we'll wait for
# star() to fill-in a value, and scale the design then        
        if (inscribed > 0)
            des = des / inscribed
    }
    
    # Add vars from left-hand side, if any
    if (inherits(basis, "formula") & (length(basis) > 2)) {
        yvars = all.vars(basis[[2]])
        for (v in yvars)  des[[v]] = NA
### not used        rsd$other=yvars
    }
    
    # Figure out sort order
###    row.names(des) = 1:nrow(des)
    des = .randomize(des, randomize=randomize)
    
    if (missing(coding)) {
        coding = sapply(xvars, function(v) as.formula(paste(v,"~",v,".as.is", sep="")))
    }
    des = as.coded.data (des, formulas=coding)
    
#     rsd$n0 = c(n0,NA)
#     rsd$non0 = c(nrow(des)-n0,NA)
#     rsd$alpha = NA
    rsd$call = CALL
#    rsd$blk.info = list(rsd[c("n0","non0","alpha")])
    attr(des, "rsdes") = rsd
    
    des
}



### new version of ccd -- constructs it using cube, star, dupe, foldover, and djoin
ccd = function(basis, generators, blocks="Block", n0=4, alpha="orthogonal", 
               wbreps=1, bbreps=1, randomize=TRUE, inscribed=FALSE, coding,
               oneblock = FALSE)
{
    n0 = rep(n0,2)[1:2]
    wbreps = rep(wbreps,2)[1:2]
    bbreps = rep(bbreps,2)[1:2]
    
# Note: We play tricks and pass a negative value if inscribed is TRUE
    if (inscribed && is.logical(inscribed)) inscribed = -999
    if (is.character(blocks)) { # need to call cube with missing blockgen
        cp = des = cube(basis, generators, n0=n0[1], reps=wbreps[1], coding=coding,
                        randomize=randomize, inscribed=inscribed)
        blkname = blocks
    }
    else if (inherits(blocks, "formula")) { # we have fractional blocks
        blkname = if(length(blocks) == 3) as.character(blocks[[2]]) else "Block"
        ### get only the rhs to make life easier for cube
        blocks = blocks[[length(blocks)]]
        cp = des = cube(basis, generators, blockgen=blocks, n0=n0[1], reps=wbreps[1], 
                        coding=coding, randomize=randomize, bid=1, inscribed=inscribed)
        # we need to add the other fractional blocks
        nfb = nrow(attr(cp, "rsdes")$blk.signs)
        if(nfb > 1) for (i in 2:nfb)
            des = djoin(des, foldover(cp, bid=i, randomize=randomize), blkname=blkname)
        cp = des
    }
    else
        stop("improper 'blocks' argument; must be a string or a formula")
    
    # Take care of replicate blocks
    if(bbreps[1] > 1) for (i in 2:bbreps[1])  {
        des = djoin(des, dupe(cp, randomize=randomize), blkname=blkname)
    }

    # Star part
    sp = star(des, n0=n0[2], alpha=alpha, reps=wbreps[2], randomize=randomize)
    # and add it the right no of times
    for (i in 1:bbreps[2])
        des = djoin(des, dupe(sp, randomize=randomize), blkname=blkname)
    
    if (inscribed) {
        div = attr(des, "rsdes")$divisor = attr(sp, "rsdes")$divisor
        pri = attr(des, "rsdes")$primary
        des[, pri] = des[, pri] / div
    }
    
    if (oneblock) {
        des[[blkname]] = NULL
        if (randomize)
            des = .randomize(des, TRUE)
    }
    
    des
}


### Make a copy of a design, re-randomizing as instructed
dupe = function(design, randomize=TRUE, coding) {
    des = .randomize(design, randomize=randomize)
    if(!missing(coding))
        codings(des) <- coding
    des
}


### New star function
# slightly older version, commented out - we're replacing with version that ac=tually 
# checks design moments
# star = function(basis, n0=4, alpha="orthogonal", reps=1, randomize=TRUE) {
#     CALL = match.call()
#     if (missing(basis)) return(CALL)
#     rsd = NULL # so these are named outside of a program block
#     if (inherits(basis, "data.frame")) {
#         rsd = attr(basis, "rsdes")
#         if (is.null(rsd))
#             stop("'basis' design must have 'rsdes' attribute.\nTry updating it via 'as.coded.data'")
#     } 
#     else
#         stop("Requires a coded.data object")
#     privars = rsd$primary
#     k = length(privars)
#     
#     # internal code to test for whether we have a cube block
#     .is.cube = function(bi) {
#         if (is.null(bi)) FALSE
#         else {
#             n0 = bi$n0
#             if (is.null(n0)) FALSE
#             else !is.na(bi$n0[1])
#         }
#     }
#     
#     if (!is.numeric(alpha)) {
#         mt = pmatch(alpha, c("orthogonal", "rotatable", "spherical", "faces"))
#         if (is.na(mt)) stop("illegal 'alpha' option")
#         got.it = FALSE
#         F.tot = 0
#         # Should be simpler than this, but this works...
#         if (is.null(rsd$blk.info)) {
#             if (.is.cube(rsd)) {
#                 f = c(rsd$non0, rsd$n0)
#                 F.tot = F.tot + rsd$non0[1]
#                 got.it = TRUE
#             }
#         }
#         else {
#             for (ann in rev(rsd$blk.info)) {
#                 if (.is.cube(ann)) {
#                     if(!got.it) f = c(ann$non0, ann$n0)
#                     F.tot = F.tot + rsd$non0[1]
#                     got.it = TRUE
#                 }
#             }
#         }
#         
#         
#         if (!got.it) stop("No cube blocks available to use in determining alpha")
#         f = f[c(1,3)]
#         alpha = switch(mt,
#                        sqrt(f[1]*(2*k*reps + n0)/(2*reps*sum(f))),
#                        sqrt(sqrt(F.tot / reps)),
#                        sqrt(k),
#                        1
#         )
#         if (!is.null(rsd$divisor)) {
#             if (rsd$divisor > 0) alpha = alpha / rsd$divisor
#             else rsd$divisor = alpha
#         }
#         
#     }
#     
#     des = as.data.frame(matrix(0, nrow=2*reps*k + n0, ncol=k))
#     names(des) = privars
#     for (i in 1:k) {
#         des[2*i-1, i] = -alpha
#         des[2*i,   i] =  alpha
#     }
#     if (reps > 1) {
#         rowrng = 1:(2*k)
#         for (i in 1:(reps-1))
#             des[rowrng + 2*k*i, ] = des[rowrng, ]
#     }
#     
#     if (is.coded.data(basis))
#         des = as.coded.data(des, formulas = codings(basis))
# 
#     des = .randomize(des, randomize=randomize)
# 
#     rsd$n0 = c(NA,n0)
#     rsd$non0 = c(NA,nrow(des)-n0)
#     rsd$alpha = alpha    
#     attr(des, "rsdes") = rsd
#     
#     des
# }


star = function(basis, n0=4, alpha="orthogonal", reps=1, randomize=TRUE) {
    CALL = match.call()
    if (missing(basis)) return(CALL)
    rsd = NULL # so these are named outside of a program block
    if (inherits(basis, "data.frame")) {
        rsd = attr(basis, "rsdes")
        if (is.null(rsd))
            stop("'basis' design must have 'rsdes' attribute.\nTry updating it via 'as.coded.data'")
    } 
    else
        stop("Requires a coded.data object")
    privars = rsd$primary
    k = length(privars)
    N = 2*k*reps + n0
    
    if (!is.numeric(alpha)) {
        mt = pmatch(alpha, c("orthogonal", "rotatable", "spherical", "faces"))
        
        if (mt == 1) { # orthogonal -- NOTE could have different alpha for each variable
            mom = .orth.moments(basis)
            if (all(is.na(mom)))
                stop("Orthogonal blocks are not achievable")
            if (any(is.na(mom))) {
                warning("Some blocks are not orthogonal first-order designs: these were ignored")
                mom = mom[!is.na(mom[,2]), , drop=FALSE]
            }
            else if (!is.matrix(mom))
                mom = matrix(mom, nrow=1)
            ratios = apply(mom[, -1, drop = FALSE], 2, function(x) x / mom[, 1])
            alpha = alphas = sqrt(N*ratios/(2*reps))
            if (is.matrix(alphas)) {
                alpha = alphas[nrow(alphas), ]
                tst = apply(alphas[-nrow(alphas), , drop=FALSE], 1, function(a) {
                    any(abs(a-alpha) > .01)
                    })
                if (any(tst))
                    warning(paste("New block is orthogonal to block", nrow(alphas),
                                  "but not to block(s)", which(tst)))
            }
        }        
        else if (mt == 2) { # rotatable
            mom = .rot.moments(basis)
            if (any(is.na(mom)))
                stop("Rotatable design is not achievable: inconsistent design moments")
            else if (is.matrix(mom))
                mom = apply(mom, 2, sum)  # overall [iiii] and [iijj] values
            alpha = (3*mom[2] - mom[1]) / (2 * reps)
            if (alpha < .0001) # after 4th root, need alpha > .1 to be valid
                stop("Rotatable design is not achievable: requires small or imaginary alpha")
            alpha = alpha^.25
        }
        else # 3 = spherical, 4 = face
            alpha = switch(mt - 2, sqrt(k), 1)
    }
        
    # possible to have several alphas -- expand so we have at least k of 'em
    if(length(alpha) < k) alpha = rep(alpha, k) 
    
    #         if (is.na(mt)) stop("illegal 'alpha' option")
    #         alpha = switch(mt,
    #                        sqrt(f[1]*(2*k*reps + n0)/(2*reps*sum(f))),
    #                        sqrt(sqrt(F.tot / reps)),
    #                        sqrt(k),
    #                        1
    #         )
    if (!is.null(rsd$divisor)) {
        if (rsd$divisor > 0) alpha = alpha / rsd$divisor
        else rsd$divisor = max(alpha)
    }
    
    des = as.data.frame(matrix(0, nrow=2*reps*k + n0, ncol=k))
    names(des) = privars
    for (i in 1:k) {
        des[2*i-1, i] = -alpha[i]
        des[2*i,   i] =  alpha[i]
    }
    if (reps > 1) {
        rowrng = 1:(2*k)
        for (i in 1:(reps-1))
            des[rowrng + 2*k*i, ] = des[rowrng, ]
    }
    
    if (is.coded.data(basis))
        des = as.coded.data(des, formulas = codings(basis))

    des = .randomize(des, randomize=randomize)
    
#     rsd$n0 = c(NA,n0)
#     rsd$non0 = c(NA,nrow(des)-n0)
#     rsd$alpha = alpha    
    rsd$call = CALL
    attr(des, "rsdes") = rsd
    
    des
}


# fold over a design. If variables is missing, use all primary ones
foldover = function(basis, variables, bid, randomize=TRUE) {
    CALL = match.call()
    if (missing(basis)) {
        return (CALL)
    }
    design = basis
    if (!inherits(design, "data.frame")) { # Try a switcheroo
        CALL$basis = NULL
        CALL$variables = basis
        return(CALL)
    }
    
    rsd = attr(basis, "rsdes")
    
    if (!missing(bid)) {
        if(is.null(rsd)) stop("information required for 'bid' is unavailable")
        signs = rsd$blk.signs
        if (is.null(signs)) stop("'bid' requires that basis be generated as a fractional block")
        s = as.numeric(signs[bid, ]) * as.numeric(signs[rsd$bid, ])
        variables = names(signs)[s < 0]
        rsd$bid = bid
    }
    else if (missing(variables)) {
        if (is.null(rsd)) 
            variables = names(design)[which(sapply(design, is.numeric))]
        else
            variables = rsd$primary
    }
    for (nm in variables)
        design[[nm]] = -design[[nm]]
    
    design = .randomize(design, randomize=randomize)
    
    attr(design, "rsdes") = rsd
    design
}

djoin = function(design1, design2, ..., blkname = "Block", blocklev) {
    if (!inherits(design1, "data.frame"))
        stop("design1 must inherit from 'data.frame'.")
    rsd1 = attr(design1, "rsdes")
    if (is.null(rsd1)) 
        stop(paste("Requires an updated 'coded.data' object\n",
            "Please update", as.character(substitute(design1)), 
                   "via 'as.coded.data' before calling"))
    
    # if design2 is actually a call, supply design1 as its 'basis' and call it
    if (is.call(design2)) {
        design2$basis = design1
        design2 = eval(design2, parent.frame())
    }
    
    # info for design 2 -- before messing with columns
    if (!is.coded.data(design2)) {
        newcode = NULL
        # if primary vars missing, make a coded dataset
        mat = match(rsd1$primary, names(design2))
        if (any(is.na(mat)))
            design2 = coded.data(design2, formulas=codings(design1))
    }
    if (is.coded.data(design2)) {
        newcode = codings(design2)
        # are they the same as design1?
        oldcode = codings(design1)
        comp = sapply(1:length(oldcode), function(i) newcode[[i]]==oldcode[[i]])
        if (!all(comp)) {
            design1 = recode.data(design1, formulas=newcode)
        }
    }
    rsd2 = attr(design2, "rsdes")
    
    # preserve annotations of separate blocks -- before modifying anything ...
#     if(is.null(rsd1$blk.info))
#         rsd1$blk.info = list(rsd1[c("n0","non0","alpha")])
#     rsd1$blk.info[[1+length(rsd1$blk.info)]] = rsd2[c("n0","non0","alpha")]

    # check names2 vs names1 - any NAs here need to be discarded in design2
    cols = match(names(design2), names(design1))
    design2 = design2[ , !is.na(cols)]
    
    # check names1 vs names2 - need to add NA columns in design2 for these
    cols = match(names(design1), names(design2))
    na.nm = names(design1)[is.na(cols)]
    for (nm in na.nm) design2[[nm]] = NA
    
# --- at this point, design1 and design2 now have the same column names
# --- (but in possibly different orders)    
    # update blocks
    blknm = rsd1$block[1]
    if (is.null(blknm) || is.na(match(blknm, names(design1)))) { # Create a blocking factor
        blknm = blkname
        design1[[blknm]] = factor(1)
        design2[[blknm]] = factor(2)
    }
    if (is.na(design2[[blknm]][1]))
        design2[[blknm]] = factor(2)
    blevs = levels(design1[[blknm]])
    nlevs = levels(design2[[blknm]])
    if (!missing(blocklev)) {
        if (length(blocklev) == length(nlevs)) {
            if (all(is.na(match(blocklev, blevs))))
                nlevs = blocklev
            else
                warning("blocklev ignored -- duplicates existing block levels")
        }
        else 
            warning("blocklev ignored -- wrong number of levels")
    }
    # redefine design2's block variable if it matches existing levels
    if (any(!is.na(match(nlevs, blevs)))) {
        # try extending a numeric sequence
        nlevs = as.character(seq(1+length(blevs), length=length(nlevs)))
        if (any(!is.na(match(nlevs, blevs))))
            nlevs = paste("A", blevs, sep="")   
    }
    levels(design2[[blknm]]) = nlevs
    

### Join the designs    
    des = rbind(design1, design2)
    
#     # update attributes
#     if (!is.null(rsd2)) {
#         for (i in 1:2) {
#             if (!is.na(rsd2$n0[i])) rsd1$n0[i] = rsd2$n0[i]
#             if (!is.na(rsd2$non0[i])) rsd1$non0[i] = rsd2$non0[i]
#         }
#         if (!is.na(rsd2$alpha)) rsd1$alpha = rsd2$alpha
#     }
    rsd1$block = blknm
    attr(des, "rsdes") = rsd1
    
    # handle dots
    for (dzn in list(...))
        des = djoin(des, dzn)
    
    des
}

stdorder = function(design) {
    if (!is.coded.data(design))
        stop("Requires a coded.data object")
    if (!("std.order" %in% names(design)))
        stop("No \"std.order\" variable found")
    if (!is.null(blknm <- attr(design, "rsdes")$block))
        key = apply(design[, c(blknm, "std.order")], 1, paste, collapse="")
    else
        key = design$std.order
    design[order(key), ]
}
