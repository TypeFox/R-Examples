### This file has summary.ref.grid S3 method and related functions

# Computes the quadratic form y'Xy after subsetting for the nonzero elements of y
.qf.non0 = function(X, y) {
    ii = (zapsmall(y) != 0)
    if (any(ii))
        sum(y[ii] * (X[ii, ii, drop = FALSE] %*% y[ii]))
    else 0
}

# utility fcn to get est's, std errors, and df
# new arg: do.se -- if FALSE, just do the estimates and return 0 for se and df
# returns a data.frame with an add'l "link" attribute if misc$tran is non-null
# .est.se.df = function(linfct, bhat, nbasis, V, dffun, dfargs, misc, do.se=TRUE, 
#                       tol=getOption("lsmeans")$estble.tol) {
# 2.13: Revised to call w/ just object instead of all those args (except linfct)
#   Also moved offset comps to here, and provided for misc$estHook
.est.se.df = function(object, do.se=TRUE, tol = get.lsm.option("estble.tol")) {
    if (nrow(object@grid) == 0) {
        result = data.frame(NA, NA, NA)
        names(result) = c(object@misc$estName, "SE", "df")
        return(result[-1, ])
    }
    misc = object@misc
    if (!is.null(hook <- misc$estHook)) {
        if (is.character(hook)) hook = get(hook)
        result = hook(object, do.se=do.se, tol=tol)
    }
    else {
        active = which(!is.na(object@bhat))
        bhat = object@bhat[active]
        result = t(apply(object@linfct, 1, function(x) {
            if (estimability::is.estble(x, object@nbasis, tol)) {
                x = x[active]
                est = sum(bhat * x)
                if(do.se) {
                    se = sqrt(.qf.non0(object@V, x))
                    df = object@dffun(x, object@dfargs)
                }
                else # if these unasked-for results are used, we're bound to get an error!
                    se = df = 0
                c(est, se, df)
            }
            else c(NA,NA,NA)
        }))
            
        if (!is.null(object@grid$.offset.))
            result[, 1] = result[, 1] + object@grid$.offset.
    }
    result = as.data.frame(result)
    names(result) = c(misc$estName, "SE", "df")

    if (!is.null(misc$tran) && (misc$tran != "none")) {
        link = if(is.character(misc$tran))
            .make.link(misc$tran)
        else if (is.list(misc$tran))
            misc$tran
        else 
            NULL
        
        if (is.list(link)) {  # See if multiple of link is requested
            if (!is.null(misc$tran.mult))
                link$mult = misc$tran.mult
            if (!is.null(link$mult))
                link = with(link, list(
                    linkinv = function(eta) linkinv(eta / mult),
                    mu.eta = function(eta) mu.eta(eta / mult) / mult,
                    name = paste0(round(mult, 3), "*", name)))
        }
        
        if (!is.null(link) && is.null(link$name))
                link$name = "linear-predictor"
        attr(result, "link") = link
    }
    result
}

# utility to compute an adjusted p value
# tail is -1, 0, 1 for left, two-sided, or right
# Note fam.info is c(famsize, ncontr, estTypeIndex)
# 2.14: added corrmat arg, dunnettx & mvt adjustments
# NOTE: corrmat is NULL unless adjust == "mvt"
.adj.p.value = function(t, df, adjust, fam.info, tail, corrmat) {
    fam.size = fam.info[1]
    n.contr = fam.info[2]
    if (n.contr == 1) # Force no adjustment when just one test
        adjust = "none"
    
    # do a pmatch of the adjust method, case insensitive
    adj.meths = c("sidak", "tukey", "scheffe", "dunnettx", "mvt", p.adjust.methods)
    k = pmatch(tolower(adjust), adj.meths)
    if(is.na(k))
        stop("Adjust method '", adjust, "' is not recognized or not valid")
    adjust = adj.meths[k]
    if ((tail != 0) && (adjust %in% c("tukey", "scheffe", "dunnettx"))) # meth not approp for 1-sided
        adjust = "sidak"
    if ((fam.info[3] != 3) && adjust == "tukey") # not pairwise
        adjust = "sidak"
    
    # asymptotic results when df is NA
    df[is.na(df)] = Inf
    
    # if estType is "prediction", use #contrasts + 1 as family size
    # (produces right Scheffe CV; Tukey ones are a bit strange)
    scheffe.dim = ifelse(fam.info[3] == 1, fam.size, fam.size - 1)
    abst = abs(t)
    if (tail == 0)
        unadj.p = 2*pt(abst, df, lower.tail=FALSE)
    else
        unadj.p = pt(t, df, lower.tail = (tail<0))
    
    if (adjust %in% p.adjust.methods) {
        if (n.contr == length(unadj.p))
            pval = p.adjust(unadj.p, adjust, n = n.contr)
        else
            pval = as.numeric(apply(matrix(unadj.p, nrow=n.contr), 2, 
                                    function(pp) p.adjust(pp, adjust, n=sum(!is.na(pp)))))
    }
    else pval = switch(adjust,
                       sidak = 1 - (1 - unadj.p)^n.contr,
                       # NOTE: tukey, scheffe, dunnettx all assumed 2-sided!
                       tukey = ptukey(sqrt(2)*abst, fam.size, zapsmall(df), lower.tail=FALSE),
                       scheffe = pf(t^2/scheffe.dim, scheffe.dim, df, lower.tail=FALSE),
                       dunnettx = 1 - .pdunnx(abst, n.contr, df),
                       mvt = 1 - .my.pmvt(t, df, corrmat, -tail) # tricky - reverse the tail because we're subtracting from 1 
                )
    chk.adj = match(adjust, c("none", "tukey", "scheffe"), nomatch = 99)
    do.msg = (chk.adj > 1) && (n.contr > 1) && 
        !((fam.size == 2) && (chk.adj < 10)) 
    if (do.msg) {
#         xtra = if(chk.adj < 10) paste("a family of", fam.size, "tests")
#         else             paste(n.contr, "tests")
        xtra = switch(adjust, 
                      tukey = paste("for comparing a family of", fam.size, "estimates"),
                      scheffe = paste("with dimensionality", scheffe.dim),
                      paste("for", n.contr, "tests")
                )
        mesg = paste("P value adjustment:", adjust, "method", xtra)
    }
    else mesg = NULL
    list(pval=pval, mesg=mesg, adjust=adjust)
}

# Code needed for an adjusted critical value
# returns a list similar to .adj.p.value
# 2.14: Added tail & corrmat args, dunnettx & mvt adjustments
# NOTE: corrmat is NULL unless adjust == "mvt"
.adj.critval = function(level, df, adjust, fam.info, tail, corrmat) {
    mesg = NULL
    
    fam.size = fam.info[1]
    n.contr = fam.info[2]
    if (n.contr == 1) # Force no adjustment when just one interval
        adjust = "none"
    
    adj.meths = c("sidak", "tukey", "scheffe", "dunnettx", "mvt", "bonferroni", "none")
    k = pmatch(tolower(adjust), adj.meths)
    if(is.na(k)) {
        k = which(adj.meths == "bonferroni") ###none")
        ###mesg = paste("\"", adjust, "\" adjustment is not valid for CIs: used bonferroni instead", sep = "")
        ###mesg = "Confidence levels are NOT adjusted for multiplicity"
    }
    adjust = adj.meths[k]
    if ((fam.info[3] != 3) && adjust == "tukey") # not pairwise
        adjust = "sidak"
    if ((tail != 0) && (adjust %in% c("tukey", "scheffe", "dunnettx"))) # meth not approp for 1-sided
        adjust = "sidak"
    if ((fam.info[3] != 3) && adjust == "tukey") # not pairwise
        adjust = "sidak"
    
    # asymptotic results when df is NA
    df[is.na(df)] = Inf
    
    scheffe.dim = ifelse(fam.info[3] == 1, fam.size, fam.size - 1)
    
    chk.adj = match(adjust, c("none", "tukey", "scheffe"), nomatch = 99)
    do.msg = (chk.adj > 1) && (n.contr > 1) && 
        !((fam.size == 2) && (chk.adj < 10)) 
    
    if (do.msg) {
#        xtra = if(chk.adj < 10) paste("a family of", fam.size, "estimates")
#        else             paste(n.contr, "estimates")
        xtra = switch(adjust, 
                      tukey = paste("for comparing a family of", fam.size, "estimates"),
                      scheffe = paste("with dimensionality", scheffe.dim),
                      paste("for", n.contr, "estimates")
        )
        mesg = paste("Conf-level adjustment:", adjust, "method", xtra)
    }
    
    adiv = ifelse(tail == 0, 2, 1) # divisor for alpha where needed
    
    cv = switch(adjust,
                none = -qt((1-level)/adiv, df),
                sidak = -qt((1 - level^(1/n.contr))/adiv, df),
                bonferroni = -qt((1-level)/n.contr/adiv, df),
                tukey = qtukey(level, fam.size, df) / sqrt(2),
                scheffe = sqrt(scheffe.dim * qf(level, scheffe.dim, df)),
                dunnettx = .qdunnx(level, n.contr, df),
                mvt = .my.qmvt(level, df, corrmat, tail)
    )
    list(cv = cv, mesg = mesg, adjust = adjust)
}


### My own functions to ease access to mvt functions
### These use one argument at a time and expands each (lower, upper) or p to a k-vector
### Use tailnum = -1, 0, or 1
### NOTE: corrmat needs "by.rows" attribute to tell which rows
###   belong to which submatrix.
.my.pmvt = function(x, df, corrmat, tailnum) {
    lower = switch(tailnum + 2, -Inf, -abs(x), x)
    upper = switch(tailnum + 2, x, abs(x), Inf)
    by.rows = attr(corrmat, "by.rows")
    if (is.null(by.rows)) 
        by.rows = list(seq_len(length(x)))
    by.sel = numeric(length(x))
    for (i in seq_along(by.rows))
        by.sel[by.rows[[i]]] = i
    df = .fix.df(df)
    apply(cbind(lower, upper, df, by.sel), 1, function(z) {
        idx = by.rows[[z[4]]]
        k = length(idx)
        pval = try(mvtnorm::pmvt(rep(z[1], k), rep(z[2], k), 
                        df = as.integer(z[3]), corr = corrmat[idx, idx]), 
                    silent = TRUE)
        if (inherits(pval, "try-error"))   NA
        else                               pval
    })
}

# Vectorized for df but needs p to be scalar
.my.qmvt = function(p, df, corrmat, tailnum) {
    tail = c("lower.tail", "both.tails", "lower.tail")[tailnum + 2] 
    df = .fix.df(df)
    by.rows = attr(corrmat, "by.rows")
    if (is.null(by.rows)) 
        by.rows = list(seq_len(length(df)))
    by.sel = numeric(length(df))
    for (i in seq_along(by.rows))
        by.sel[by.rows[[i]]] = i
    # If df all equal, compute just once for each by group
    eq.df = (diff(range(df)) == 0)
    i1 = if (eq.df)   sapply(by.rows, function(r) r[1])
         else         seq_along(df)
    result = apply(cbind(p, df[i1], by.sel[i1]), 1, function(z) {
        idx = by.rows[[z[3]]]
        cv = try(mvtnorm::qmvt(z[1], tail = tail, 
                    df = as.integer(z[2]), corr = corrmat[idx, idx])$quantile,
                 silent = TRUE)
        if (inherits(cv, "try-error"))     NA
        else                               cv
    })
    if (eq.df) {
        res = result
        result = numeric(length(df))
        for(i in seq_along(by.rows))
            result[by.rows[[i]]] = res[i]
    } 
    result
}

# utility to get appropriate integer df
.fix.df = function(df) {
    sapply(df, function(d) {
        if (d > 0) d = max(1, d)
        if (is.infinite(d) || (d > 9999)) d = 0
        floor(d + .25) # tends to round down
    })
}

### My approximate dunnett distribution 
### - a mix of the Tukey cdf and Sidak-corrected t
.pdunnx = function(x, k, df, twt = (k - 1)/k) {
    tukey = ptukey(sqrt(2)*x, (1 + sqrt(1 + 8*k))/2, df)
    sidak = (pf(x^2, 1, df))^k
    twt*tukey + (1 - twt)*sidak
}

# Uses linear interpolation to get quantile
.qdunnx = function(p, k, df, ...) {
     if (k < 1.005)
         return(qt(1 - .5*(1 - p), df))
    xtuk = qtukey(p, (1 + sqrt(1 + 8*k))/2, df) / sqrt(2)
    xsid = sqrt(qf(p^(1/k), 1, df))
    fcn = function(x, d) 
        .pdunnx(x, k, d, ...) - p
    apply(cbind(xtuk, xsid, df), 1, function(r) {
        if (abs(diff(r[1:2])) < .0005)
            return (r[1])
        x = try(uniroot(fcn, r[1:2], tol = .0005, d = r[3]), silent = TRUE)
        if (inherits(x, "try-error")) {
            warning("Root-finding failed; using qtukey approximation for Dunnett quantile")
            return(xtuk)
        }
        else
            x$root
    })
}



### Support for different prediction types ###

# Valid values for type arg or predict.type option
.valid.types = c("link","response","lp","linear")

# get "predict.type" option from misc, and make sure it's legal
.get.predict.type = function(misc) {
    type = misc$predict.type
    if (is.null(type))
        .valid.types[1]
    else
        .validate.type(type)
}

# check a "type" arg to make it legal
.validate.type = function (type) {
    .valid.types[pmatch(type, .valid.types, 1)]
}

# S3 predict method
predict.ref.grid <- function(object, type, ...) {
    # update with any "summary" options
    opt = get.lsm.option("summary")
    if(!is.null(opt)) {
        opt$object = object
        object = do.call("update.ref.grid", opt)
    }
    
    if (missing(type))
        type = .get.predict.type(object@misc)
    else
        type = .validate.type(type)
    
    pred = .est.se.df(object, do.se=FALSE)
    result = pred[[1]]
    # MOVED TO .EST.SE.DF    
    #     if (".offset." %in% names(object@grid))
    #         result = result + object@grid[[".offset."]]
    if (type == "response") {
        link = attr(pred, "link")
        if (!is.null(link))
            result = link$linkinv(result)
    }
    result
}

# S3 summary method
summary.ref.grid <- function(object, infer, level, adjust, by, type, df, 
                             null, delta, side, ...) {
    ### For missing arguments, get from misc, else default    
    if(missing(infer))
        infer = object@misc$infer
    if(missing(level))
        level = object@misc$level
    if(missing(adjust))
        adjust = object@misc$adjust
    if(missing(by))
        by = object@misc$by.vars
    
    if (missing(type))
        type = .get.predict.type(object@misc)
    else
        type = .validate.type(type)
    
    if(missing(df)) 
        df = object@misc$df
    if(!is.null(df))
        object@dffun = function(k, dfargs) df
    
    # for missing args that default to zero unless provided or in misc slot
    .nul.eq.zero = function(val) {
        if(is.null(val)) 0
        else val
    }
    if(missing(null))
        null = .nul.eq.zero(object@misc$null)
    if(missing(delta))
        delta = .nul.eq.zero(object@misc$delta)
    if(missing(side))
        side = .nul.eq.zero(object@misc$side)
    
    # update with any "summary" options
    opt = get.lsm.option("summary")
    if(!is.null(opt)) {
        opt$object = object
        object = do.call("update.ref.grid", opt)
    }
    
    
    # reconcile all the different ways we could specify the alternative
    # ... and map each to one of the first 3 subscripts
    side.opts = c("left","both","right","two-sided","noninferiority","nonsuperiority","equivalence","superiority","inferiority","0","2","-1","1","+1","<",">","!=","=")
    side.map =  c( 1,     2,     3,      2,          3,               1,               2,            3,            1,            2,  2,   1,  3,   3,  1,  3,  2,   2)
    side = side.map[pmatch(side, side.opts, 2)[1]] - 2
    delta = abs(delta)
    
    result = .est.se.df(object)
    
    lblnms = setdiff(names(object@grid), 
                     c(object@roles$responses, ".offset.", ".wgt."))
    lbls = object@grid[lblnms]
    
    zFlag = (all(is.na(result$df)))
    inv = (type == "response") # flag to inverse-transform
    link = attr(result, "link")
    if (inv && is.null(link))
        inv = FALSE

    if ((length(infer) == 0) || !is.logical(infer)) 
        infer = c(FALSE, FALSE)
    if(length(infer == 1)) 
        infer = c(infer,infer)
    
    if(inv && !is.null(object@misc$tran)) {
        if (!is.null(object@misc$inv.lbl))
            names(result)[1] = object@misc$inv.lbl
        else
            names(result)[1] = "lsresponse"
    }

    attr(result, "link") = NULL
    estName = names(result)[1]
    
    mesg = object@misc$initMesg
    
    ### Add an annotation when we show results on lp scale and
    ### there is a transformation
    if (!inv && !is.null(link)) {
        mesg = c(mesg, paste("Results are given on the", link$name, "(not the response) scale."))
    }
    if (inv && !is.null(link$unknown)) {
        mesg = c(mesg, paste0('Unknown transformation "', link$name, '": no transformation done'))
        inv = FALSE
        link = NULL
    }
    
    # et = 1 if a prediction, 2 if a contrast (or unmatched or NULL), 3 if pairs
    et = pmatch(c(object@misc$estType, "c"), c("prediction", "contrast", "pairs"), nomatch = 2)[1]
    
    by.size = nrow(object@grid)
    if (!is.null(by))
        for (nm in by)
            by.size = by.size / length(unique(object@levels[[nm]]))
    fam.info = c(object@misc$famSize, by.size, et)
    cnm = NULL
    
    # get vcov matrix only if needed (adjust == "mvt")
    corrmat = NULL
    if (!is.na(pmatch(adjust, "mvt"))) {
        corrmat = cov2cor(vcov(object))
        attr(corrmat, "by.rows") = .find.by.rows(object@grid, by)
    }
    
    if(infer[1]) { # add CIs
        acv = .adj.critval(level, result$df, adjust, fam.info, side, corrmat)
        ###adjust = acv$adjust # in older versions, I forced same adj method for tests
        cv = acv$cv
        cv = switch(side + 2, cbind(-Inf, cv), cbind(-cv, cv), cbind(-cv, Inf))
        cnm = if (zFlag) c("asymp.LCL", "asymp.UCL") else c("lower.CL","upper.CL")
        result[[cnm[1]]] = result[[1]] + cv[, 1]*result$SE
        result[[cnm[2]]] = result[[1]] + cv[, 2]*result$SE
        mesg = c(mesg, paste("Confidence level used:", level), acv$mesg)
        if (inv) {
            clims = with(link, cbind(linkinv(result[[cnm[1]]]), linkinv(result[[cnm[2]]])))
            idx = if (all(clims[ ,1] <= clims[, 2])) 1:2 else 2:1
            result[[cnm[1]]] = clims[, idx[1]]
            result[[cnm[2]]] = clims[, idx[2]]
            mesg = c(mesg, paste("Intervals are back-transformed from the", link$name, "scale"))
        }
    }
    if(infer[2]) { # add tests
        if (!all(null == 0)) {
            result[["null"]] = null
            if (!is.null(link))
                result[["null"]] = link$linkinv(result[["null"]])
        }
        tnm = ifelse (zFlag, "z.ratio", "t.ratio")
        tail = ifelse(side == 0, -sign(abs(delta)), side)
        if (side == 0) {
            if (delta == 0) # two-sided sig test
                t.ratio = result[[tnm]] = (result[[1]] - null) / result$SE
            else
                t.ratio = result[[tnm]] = (abs(result[[1]] - null) - delta) / result$SE
        }
        else {
            t.ratio = result[[tnm]] = (result[[1]] - null + side * delta) / result$SE            
        }
        apv = .adj.p.value(t.ratio, result$df, adjust, fam.info, tail, corrmat)
        adjust = apv$adjust   # in case it was abbreviated
        result$p.value = apv$pval
        mesg = c(mesg, apv$mesg)
        if (delta > 0)
            mesg = c(mesg, paste("Statistics are tests of", c("nonsuperiority","equivalence","noninferiority")[side+2],
                                 "with a threshold of", delta))
        if(tail != 0) 
            mesg = c(mesg, paste("P values are ", ifelse(tail<0,"left-","right-"),"tailed", sep=""))
        if (!is.null(link)) 
            mesg = c(mesg, paste("Tests are performed on the", link$name, "scale"))
    }
    if (inv) {
        result[["SE"]] = with(link, abs(mu.eta(result[[1]]) * result[["SE"]]))
        result[[1]] = with(link, linkinv(result[[1]]))
    }
    
    if (length(object@misc$avgd.over) > 0) {
        qual = attr(object@misc$avgd.over, "qualifier")
        if (is.null(qual)) qual = ""
        mesg = c(paste0("Results are averaged over", qual, " the levels of: ",
                       paste(object@misc$avgd.over, collapse = ", ")), mesg)
    }
    summ = cbind(lbls, result)
    attr(summ, "estName") = estName
    attr(summ, "clNames") = cnm  # will be NULL if infer[1] is FALSE
    attr(summ, "pri.vars") = setdiff(union(object@misc$pri.vars, object@misc$by.vars), by)
    attr(summ, "by.vars") = by
    attr(summ, "mesg") = unique(mesg)
    class(summ) = c("summary.ref.grid", "data.frame")
    summ
}


# left-or right-justify column labels for m depending on "l" or "R" in just
.just.labs = function(m, just) {
    nm = dimnames(m)[[2]]
    for (j in seq_len(length(nm))) {
        if(just[nm[j]] == "L") 
            nm[j] = format(nm[j], width = nchar(m[1,j]), just="left")
    }
    dimnames(m) = list(rep("", nrow(m)), nm)
    m
}

# Format a data.frame produced by summary.ref.grid
print.summary.ref.grid = function(x, ..., digits=NULL, quote=FALSE, right=TRUE) {
    x.save = x
    if (!is.null(x$df)) x$df = round(x$df, 2)
    if (!is.null(x$t.ratio)) x$t.ratio = round(x$t.ratio, 3)
    if (!is.null(x$z.ratio)) x$z.ratio = round(x$z.ratio, 3)
    if (!is.null(x$p.value)) {
        fp = x$p.value = format(round(x$p.value,4), nsmall=4, sci=FALSE)
        x$p.value[fp=="0.0000"] = "<.0001"
    }
    just = sapply(x.save, function(col) if(is.numeric(col)) "R" else "L")
    xc = as.matrix(format.data.frame(x, digits=digits, na.encode=FALSE))
    m = apply(rbind(just, names(x), xc), 2, function(x) {
        w = max(sapply(x, nchar))
        if (x[1] == "R") format(x[-seq_len(2)], width = w, justify="right")
        else format(x[-seq_len(2)], width = w, justify="left")
    })
    if(!is.matrix(m)) m = t(as.matrix(m))
    by.vars = attr(x, "by.vars")
    if (is.null(by.vars)) {
        m = .just.labs(m, just)
        print(m, quote=FALSE, right=TRUE)
        cat("\n")
    }
    else { # separate listing for each by variable
        m = .just.labs(m[, setdiff(names(x), by.vars)], just)
        pargs = as.list(x[,by.vars, drop=FALSE])
        pargs$sep = ", "
        lbls = do.call(paste, pargs)
        for (lb in unique(lbls)) {
            rows = which(lbls==lb)
            levs = paste(by.vars, "=", xc[rows[1], by.vars])
            cat(paste(paste(levs, collapse=", ")), ":\n", sep="")
            print(m[rows, , drop=FALSE], ..., quote=quote, right=right)
            cat("\n")
        }
    }
    
    msg = unique(attr(x, "mesg"))
    if (!is.null(msg))
        for (j in seq_len(length(msg))) cat(paste(msg[j], "\n"))
    
    invisible(x.save)
}

