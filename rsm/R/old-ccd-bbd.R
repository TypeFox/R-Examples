### This source file contains the versions of ccd and bbd from rsm version 1.41
### They are provided to allow for working around differences in rsm starting 
### with 2.00, and to allow one to exactly reproduce a design that was created
### in an old version of rsm -- including the randomization if the same seed
### is used. The code here is NOT quite identical to the old code, due to
### changes in coded.data structure. We even reproduce the old coded structure.
### These functions are NOT exported.

### generate a CCD
# basis - formula - lhs (opt): dep var name(s);  rhs: var names for basic grid
# 

.ccd.1.41 = function(basis, generators, blocks="Block", n0=4, alpha="orthogonal", 
               wbreps=1, bbreps=1, randomize=TRUE, inscribed=FALSE, coding,
               new.style = FALSE)
{
    if (inherits(basis, "formula"))
        xvars = all.vars(basis[[length(basis)]])
    else if (is.numeric(basis))
        xvars = paste("x", 1:basis, sep="")
    else
        stop("'basis' must be an integer or a formula")
    
    args = lapply(xvars, function(nm) c(-1,1))
    names(args) = xvars
    cube = do.call(expand.grid, args)
    
    if (!missing(generators)) {
        if (!is.list(generators)) generators = list(generators)
        for (gen in generators) {
            gen = as.character(gen)        
            cube[[gen[[2]]]] = with(cube, eval (parse (text = as.character(gen[[3]]))))
        }
    }
    
    k = ncol(cube)
    
    ### At first, star will be face-centered (as if alpha = 1); we'll scale it later
    star = as.data.frame (matrix(0, nrow = 2*k, ncol = k))
    xvars = names(star) = names (cube)
    for (j in 1:k) star[c(2*j-1,2*j), j] = c(-1, 1)
    
    if (length(wbreps) == 1) wbreps = rep(wbreps, 2)
    if (length(bbreps) == 1) bbreps = rep(bbreps, 2)
    
    # within-block reps
    if (wbreps[1] > 1) cube = cube[rep(1:nrow(cube), wbreps[1]), ]
    if (wbreps[2] > 1) star = star[rep(1:nrow(star), wbreps[2]), ]
    
    # Fractional blocking
    if (is.character(blocks)) {
        blknm = blocks
        nblev = 1
        blk = rep(1, nrow(cube))
        chkterm = ""
    }
    else if (inherits(blocks, "formula")) {
        blknm = as.character(blocks[[2]])
        what = as.character(blocks[[3]][[1]])
        if (what == "*") 
            gens = as.character(blocks[3])
        else 
            gens = as.character(blocks[[3]])[-1]
        bgen = lapply(gens, function(g) with(cube, eval(parse(text=g))))
        blk = as.numeric(factor(do.call(paste, bgen)))
        nblev = max(blk)
        chkterm = "factor(blk) + "
    }
    else
        stop("'blocks' must be a string or a formula")
    
    # Check for aliasing
    v = paste(names(cube), collapse=",")
    fake.resp = rnorm(nrow(cube))
    fstg = paste("fake.resp ~", chkterm, "FO(", v, ") + TWI(", v, ")")
    modl = lm(formula(fstg), data=cube)
    if (any(is.na(coef(modl))))
        warning("Some 1st or 2nd-order terms are aliased in the cube portion of this design")
    
    
    # center points
    zero = as.data.frame(matrix(rep(0, k), nrow=1))
    names(zero) = names(cube)
    if (length(n0) == 1) n0 = c(n0, n0)
    if (n0[1] > 0) {
        cube = rbind(cube, zero[rep(1, nblev*n0[1]), ])
        blk = c(blk, rep(unique(blk), n0[1]))
    }
    if (n0[2] > 0) 
        star = rbind(star, zero[rep(1, n0[2]), ])
    
    # Block reps
    nc = nrow(cube)
    if (bbreps[1] > 1) {
        cube = cube[rep(1:nc, bbreps[1]), ]
        blk = nblev*rep(0:(bbreps[1]-1), rep(nc, bbreps[1])) 
        + rep(blk, bbreps[1])
        nblev = max(blk)
    }
    ns = nrow(star)
    if (bbreps[2] > 1) star = star[rep(1:ns, bbreps[2]), ]
    sblk = rep( (1+nblev):(bbreps[2]+nblev), rep(ns, bbreps[2]))
    
    # Figure out alpha if given as criterion
    if (is.character(alpha)) {
        c.ii = sum(cube[[1]]^2)
        s.ii = sum(star[[1]]^2)
        what = pmatch(alpha, c("rotatable", "orthogonal"))
        if (is.na(what)) stop ("alpha must be 'rotatable', 'orthogonal', or a value")
        if (what==1)
            alpha = (2 * c.ii / s.ii) ^ .25
        else
            alpha = sqrt(nrow(star) / s.ii  * c.ii / nrow(cube))
    }
    if (inscribed)
        cube = cube / alpha
    else
        star = star * alpha
    
    # append blocking variable
    cube = cbind(blk, cube)
    star = cbind(sblk, star)
    names(cube)[1] = names(star)[1] = blknm
    
    # Figure out row names
    ord = order(blk, 1:nrow(cube))
    cube = cube[ord, ]
    blk = blk[ord]
    row.names(cube) = paste("C", blk, ".", rep(1:(nrow(cube)/nblev), nblev), sep="")
    row.names(star) = paste("S", sblk, ".", rep(1:(nrow(star)/bbreps[2]), bbreps[2]), sep="")
    
    # assemble design
    des = rbind(cube, star)
    
    # Add vars from left-hand side, if any
    if (inherits(basis, "formula") & (length(basis) > 2)) {
        yvars = all.vars(basis[[2]])
        for (v in yvars)  des[[v]] = NA
    }
    
    
    stdord = 1:nrow(des)
    
    # Figure out sort order
    if (randomize) {
        ord = order(des[[1]] + runif(nrow(des)))
        des = des[ord, ]
        stdord = stdord[ord]
    }
    
    des[[1]] = factor(des[[1]])
    
    if (!missing(coding)) {
        des = as.coded.data (des, formulas=coding)
        if (!new.style)
            attr(des, "rsdes") = NULL # remove rsm2.00-style attributes
    } 
    
    if (new.style) {
        if (missing(coding))
            coding = sapply(xvars, function(v) as.formula(paste(v,"~",v,".as.is", sep="")))
        des = .randomize(as.coded.data(des, formulas=coding), randomize=FALSE)
        des$std.order = stdord
    }
    
    des
}



### Generate a Box-Behnken design (up thru k=7)

.bbd.1.41 = function(k, n0=4, block = (k==4|k==5), randomize=TRUE, coding, 
                     new.style=FALSE) 
{
    args = list(fcn=bbd, k=k, n0=n0, block=block, 
                randomize=randomize, strip = !new.style)
    if (!missing(coding)) args$coding = coding
    des = do.call(.gen.des.old.way, args)
    
    if(!new.style) {
        attr(des, "rsdes") = NULL
        if (missing(coding)) {
            attr(des, "codings") = NULL
            class(des) = setdiff(class(des), "coded.data")
        }
    }
    
    des
}

### Call a new design-generating function but randomize it the old way
### Returns new-style coded data if fcn does
### strip determines whether to strip std.order and run.order
.gen.des.old.way = function(fcn, ..., randomize=TRUE, strip=FALSE) {
    des = fcn(..., randomize=FALSE)
    if (is.null(attr(des, "rsdes"))) # Not a new-style result
        return(des)
    if (randomize) {
        blknm = attr(des, "rsdes")$block
        if (is.null(blknm))  blks = 0
        else                 blks = as.numeric(des[[blknm]])
        ord = order(blks + runif(nrow(des)))
        des = des[ord, ]
        if (!strip) 
            row.names(des) = 1:nrow(des)
        ub = unique(blks)
        if (length(ub) == 1)
            des$run.order = 1:nrow(des)
        for (b in 1:length(ub))
            des$run.order[blks==b] = 1:sum(blks==b)
    }
    if (strip) {
        xcols = match(c("run.order","std.order"), names(des))
        if(length(xcols) > 0) des = des[, -xcols]
    }
    des
}
