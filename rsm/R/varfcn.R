# Variance function of a design

# Note on default vectors: if model contains call to FO, SO, PQ, TWI, or poly,
# we use variables listed in 1st such call. Else we use all non-factor predictors

varfcn = function(design, formula, dist=seq(0,2,by=.1), 
                  vectors, contour=FALSE, plot=TRUE, main, ...) 
{
    tt = delete.response(terms(formula))
    mf = model.frame(tt, design)
    mm = model.matrix(tt, mf)
    nxpxinv = nrow(mm) * solve(t(mm) %*% mm)
    
    if (missing(vectors)) {
        dc = attr(attr(mf, "terms"), "dataClasses")
        nm = names(dc)
        rsmterms = grep("FO\\(|SO\\(|TWI\\(|PQ\\(|poly\\(", nm)
        if (length(rsmterms) > 0) {
            nm = sapply(as.formula(paste("~",nm[rsmterms[1]]))[[2]], as.character)[-1]
        }
        else {
            facs = grep("factor|ordered", dc)
            if (length(facs)>0) 
                nm = nm[-facs]
        }
        vectors = as.data.frame(matrix(1, nrow=length(nm), ncol=length(nm)))
        names(vectors) = nm
        if (length(nm)>1) for (i in 2:length(nm)) vectors[[i]] [1:(i-1)] = 0
    }
    else { 
        # make sure there are no zero rows
        nz = apply(vectors, 1, function(x) sum(abs(x)>1e-4) > 0)
        vectors = vectors[nz, ]
    }
    
    if (missing(main))
        main = paste(paste(as.character(substitute(design)),collapse=""), paste(formula, collapse=" "), sep=": ")
    
    if (contour) {
        temp = sort(c(-dist[-1], dist))
        X = expand.grid(temp, temp)
        names(X) = names(vectors)[1:2]
    }
    else {
        if (ncol(vectors) > 1) 
            temp = apply(vectors, 1, function(row) row / sqrt(sum(row^2)))
        else temp=t(as.matrix(vectors))
        X = apply(temp, 1, function(vec) sapply(vec, function(v) v*dist))
        X = as.data.frame(X)
    }
    
    # Typical observation of a given variable
    typical = function(var) {
        if (is.numeric(var)) mean(var, na.rm=TRUE) else var[1]
    }
    
    # assemble frame for predictions
    n = nrow(X)
    pd = lapply(names(design), function(nm) 
        if (is.null(X[[nm]])) rep(typical(design[[nm]]), n) else X[[nm]])
    pd = as.data.frame(pd)
    names(pd) = names(design)
    pf = model.frame(tt, pd)
    pm = model.matrix(tt, pf)
    X$VF = apply(pm, 1, function(x) sum(x * (nxpxinv %*% x)))
    if (contour && plot) {
        z = matrix(X$VF, nrow=length(temp))
        nm = names(X)[1:2]
        contour(temp, temp, z, xlab=nm[1], ylab=nm[2], main=main, ...)
        tbl = with(design, do.call("table", list(get(nm[1]),get(nm[2]))))
        tmp = do.call("expand.grid", lapply(dimnames(tbl), as.numeric))
        tmp$freq = as.numeric(tbl)
        tmp = tmp[tmp$freq > 0, ]
        points(tmp[[1]], tmp[[2]], pch=18, cex=sqrt(tmp$freq), col="red")
    }
    else {
        X = cbind(dir=rep(1:nrow(vectors), each=length(dist)), dist=rep(dist, nrow(vectors)), X)
        if(plot) {
            y = matrix(X$VF, nrow=length(dist))
            mpargs = list(x=dist, y=y, type="l", xlab="Distance from center", 
                          ylab="Scaled prediction variance", main=main, lty=1)
            # Let user override my choice of lty
            if (!is.null(list(...)$lty)) 
                mpargs$lty = list(...)$lty
            do.call(matplot, mpargs)
        }
    }
    invisible(X)
}
