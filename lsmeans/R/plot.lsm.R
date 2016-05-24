# S3 plot method for lsmobj objects (NOT ref.grid as relies on pri.fac attribute etc.)
# ... are arguments sent to update()


plot.lsmobj = function(x, y, type, intervals = TRUE, comparisons = FALSE, 
                       alpha = .05, adjust = "tukey", int.adjust = "none", ...) {
    if(!missing(type))
        object = update(x, predict.type = type, ..., silent = TRUE)
    else
        object = update(x, ..., silent = TRUE)
    if (missing(int.adjust)) {
        int.adjust = object@misc$adjust
        if (is.null(int.adjust))
            int.adjust = "none"
    }
    summ = summary(object, infer = c(TRUE, FALSE), adjust = int.adjust)
    estName = attr(summ, "estName")
    extra = NULL
    if(comparisons) {
        extra = object
        extra@misc$comp.alpha = alpha
        extra@misc$comp.adjust = adjust
    }
    .plot.srg(x=summ, intervals = intervals, extra = extra, ...)
}

# May use in place of plot.lsmobj but no control over level etc.
# extra is a placeholder for comparison-interval stuff
plot.summary.ref.grid = function(x, y, horizontal = TRUE, xlab, ylab, layout, ...) {
    .plot.srg (x, y, horizontal, xlab, ylab, layout, ...)
}

# Workhorse for plot.summary.ref.grid
.plot.srg = function(x, y, horizontal = TRUE, xlab, ylab, layout, intervals = TRUE, extra = NULL, ...) {
        
    if (!requireNamespace("lattice"))
        stop("This function requires the 'lattice' package be installed.")
    
    summ = x # so I don't get confused
    estName = attr(summ, "estName")
    clNames = attr(summ, "clNames")
    if (is.null(clNames)) {
        warning("No information available to display confidence limits")
        lcl = ucl = summ[[estName]]
    }
    else {
        lcl = summ[[clNames[1]]]
        ucl = summ[[clNames[2]]]
    }
    
    # Panel functions...
    prepanel.ci = function(x, y, horizontal=TRUE, intervals=TRUE,
                           lcl, ucl, subscripts, ...) {
        x = as.numeric(x)
        lcl = as.numeric(lcl[subscripts])
        ucl = as.numeric(ucl[subscripts])
        if (!intervals) # no special scaling needed
            list()
        else if (horizontal)
            list(xlim = range(x, ucl, lcl, finite = TRUE)) 
        else
            list(ylim = range(y, ucl, lcl, finite = TRUE)) 
    }
    panel.ci <- function(x, y, horizontal=TRUE, intervals=TRUE,
                         lcl, ucl, lcmpl, rcmpl,                          subscripts, pch = 16, 
                         lty = dot.line$lty, lwd = dot.line$lwd, 
                         col = dot.symbol$col, col.line = dot.line$col, ...) {
        dot.line <- lattice::trellis.par.get("dot.line")
        dot.symbol <- lattice::trellis.par.get("dot.symbol")
        x = as.numeric(x)
        y = as.numeric(y)
        lcl = as.numeric(lcl[subscripts])
        ucl = as.numeric(ucl[subscripts])
        compare = !is.null(lcmpl)
        if(compare) {
            lcmpl = as.numeric(lcmpl[subscripts])
            rcmpl = as.numeric(rcmpl[subscripts])
        }
        if(horizontal) {
            lattice::panel.abline(h = unique(y), col = col.line, lty = lty, lwd = lwd)
            if(intervals) 
                lattice::panel.arrows(lcl, y, ucl, y, col = col, length = .6, unit = "char", angle = 90, code = 3)
            if(compare) {
                s = (x > min(x))
                lattice::panel.arrows(lcmpl[s], y[s], x[s], y[s], length = .5, unit = "char", code = 1, col = "red", type = "closed", fill="red")
                s = (x < max(x))
                lattice::panel.arrows(rcmpl[s], y[s], x[s], y[s], length = .5, unit = "char", code = 1, col = "red", type = "closed", fill="red")
            }
        }
        else {
            lattice::panel.abline(v = unique(x), col = col.line, lty = lty, lwd = lwd)
            if(intervals)
                lattice::panel.arrows(x, lcl, x, ucl, col=col, length = .6, unit = "char", angle = 90, code = 3)
            if(compare) {
                s = (y > min(y))
                lattice::panel.arrows(x[s], lcmpl[s], x[s], y[s], length = .5, unit = "char", code = 1, col = "red", type = "closed", fill="red")
                s = (y < max(y))
                lattice::panel.arrows(x[s], rcmpl[s], x[s], y[s], length = .5, unit = "char", code = 1, col = "red", type = "closed", fill="red")
            }
        }
        lattice::panel.xyplot(x, y, pch=16, ...)
    }
    my.strip = lattice::strip.custom(strip.names = c(TRUE,TRUE), strip.levels = c(TRUE,TRUE), sep = " = ")
    
    priv = attr(summ, "pri.vars")
    pf = do.call(paste, summ[priv])
    summ$pri.fac = factor(pf, levels=unique(pf))
    chform = ifelse(horizontal,
                    paste("pri.fac ~", estName),
                    paste(estName, "~ pri.fac"))
    
    byv = attr(summ, "by.vars")
    if (!is.null(byv)) {
        chform = paste(chform, "|", paste(byv, collapse="*"))
        lbv = do.call("paste", summ[byv]) # strings for matching by variables
        ubv = unique(lbv)
    }
    else {
        lbv = rep(1, nrow(summ))
        ubv = 1
    }
    
    
    # Obtain comparison limits
    if (!is.null(extra)) {
        # we need to work on the linear predictor scale
        # typeid = 1 -> response, 2 -> other
        typeid = pmatch(extra@misc$predict.type, "response", nomatch = 2)
        if(length(typeid) < 1) typeid = 2        
        if (typeid == 1)
            est = predict(extra, type = "lp")
        else
            est = summ[[estName]]
        
        alpha = extra@misc$comp.alpha
        adjust = extra@misc$comp.adjust
        psumm = confint(pairs(extra), level = 1 - alpha, type = "lp", adjust = adjust)
        k = ncol(psumm)
        del = (psumm[[k]] - psumm[[k-1]]) / 4 # half the halfwidth, on lp scale
        diff = psumm[[attr(psumm, "estName")]]
        overlap = apply(psumm[ ,(k-1):k], 1, function(x) 2*min(-x[1],x[2])/(x[2]-x[1]))
        
        # figure out by variables and indexes (lbv, ubv already defined)
        if(is.null(byv))
            pbv = rep(1, nrow(psumm))
        else
            pbv = do.call("paste", psumm[byv])
        neach = length(lbv) / length(ubv)
        # indexes for pairs results -- est[id1] - est[id2]
        id1 = rep(seq_len(neach-1), rev(seq_len(neach-1)))
        id2 = unlist(sapply(seq_len(neach-1), function(x) x + seq_len(neach-x)))
        # list of psumm row numbers involved in each summ row
        involved = lapply(seq_len(neach), function(x) union(which(id2==x), which(id1==x)))
        
        # initialize arrays
        mind = numeric(length(lbv))   # for minima of del
        llen = rlen = numeric(neach)  # for left and right arrow lengths
        npairs = length(id1)
        iden = diag(rep(1, 2*neach))
        
        for (by in ubv) {
            d = del[pbv == by]
            rows = which(lbv == by)
            for(i in seq_len(neach)) 
                mind[rows[i]] = min(d[involved[[i]]])
            
            # Set up regression equations to match arrow overlaps with interval overlaps
            # We'll add rows later (with weights 1) to match with mind values
            lmat = rmat = matrix(0, nrow = npairs, ncol = neach)
            y = numeric(npairs)
            v1 = 1 - overlap[pbv == by]
            dif = diff[pbv == by]
            for (i in seq_len(npairs)) {
                #wgt = 6 * max(0, ifelse(v1[i] < 1, v1[i], 2-v1[i]))
                wgt = 3 + 20 * max(0, .5 - (1 - v1[i])^2)
                # really this is sqrt of weight
                if (dif[i] > 0)   # id2  <----->  id1
                    lmat[i, id1[i]] = rmat[i, id2[i]] = wgt*v1[i]
                else  # id1  <----->  id2
                    rmat[i, id1[i]] = lmat[i, id2[i]] = wgt*v1[i]
                y[i] = wgt * abs(dif[i])
            }
            X = rbind(cbind(lmat, rmat),iden)
            y = c(y, rep(mind[rows], 2))
            soln = qr.coef(qr(X), y)
            ll = llen[rows] = soln[seq_len(neach)]
            rl = rlen[rows] = soln[neach + seq_len(neach)]
            
            # Perhaps put some kind of a check here?
            for (i in seq_len(npairs)) {
                v = 1 - v1[i]
                obsv = 1 - abs(dif[i]) / ifelse(dif[i] > 0, 
                                ll[id1[i]] + rl[id2[i]], 
                                rl[id1[i]] + ll[id2[i]])
                if (v*obsv < 0)
                    message("Comparison discrepancy in group ", by, 
                            ", ", psumm[i, 1], 
                            ":\n    Target overlap = ", round(v, 4),
                            ", overlap on graph = ", round(obsv, 4))
            }
        }
        # shorten arrows that go past the data range
        rng = range(est)
        ii = which(est - llen < rng[1])
        llen[ii] = est[ii] - rng[1]
        ii = which(est + rlen > rng[2])
        rlen[ii] = rng[2] - est[ii]
        
        invtran = I
        if (typeid == 1) {
            tran = extra@misc$tran
            if(is.character(tran)) {
                link = try(make.link(tran), silent=TRUE)
                if (!inherits(link, "try-error"))
                    invtran = link$linkinv
            }
            else if (is.list(tran))
                invtran = tran$linkinv
        }
        
        lcmpl = invtran(est - llen)
        rcmpl = invtran(est + rlen)
    }
    else lcmpl = rcmpl = NULL
    
    
    if (missing(layout)) {
        layout = c(1, length(ubv))
        if(!horizontal) 
            layout = rev(layout)
    }
    
    facName = paste(priv, collapse=":")
    form = as.formula(chform)
    if (horizontal) {
        if (missing(xlab)) xlab = estName
        if (missing(ylab)) ylab = facName
        lattice::dotplot(form, prepanel=prepanel.ci, panel=panel.ci, 
                strip = my.strip, horizontal = TRUE,
                ylab = ylab, xlab = xlab,
                data = summ, intervals = intervals, lcl=lcl, ucl=ucl, 
                lcmpl=lcmpl, rcmpl=rcmpl, layout = layout, ...)
    }
    else {
        if (missing(xlab)) xlab = facName
        if (missing(ylab)) ylab = estName
        lattice::dotplot(form, prepanel=prepanel.ci, panel=panel.ci, 
                strip = my.strip, horizontal = FALSE,
                xlab = xlab, ylab = ylab,
                data = summ, intervals = intervals, lcl=lcl, ucl=ucl, 
                lcmpl=lcmpl, rcmpl=rcmpl, layout = layout, ...)
    }
}
