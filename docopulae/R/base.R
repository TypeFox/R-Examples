

#' Sequence Generation
#'
#' \code{seq1} is similar to \code{base::seq}, however \code{by} is strictly \code{1} by default and \code{integer(0)} is returned if range is empty.
#'
#' @param from,to,by see \code{\link[base]{seq}}.
#'
#' @return \code{seq1} returns either \code{integer(0)} if range is empty or what an appropriate call to \code{base::seq} returns otherwise.
#'
#' See examples below.
#'
#' @example R/examples/seq1.R
#'
#' @seealso \code{\link[base]{seq}}
#'
#' @export
seq1 = function(from, to, by=1) {
    if (to == from)
        return(from)
    if (sign(to - from) != sign(by))
        return(integer(0))
    return(seq(from, to, by))
}


Derivf = function(f, names) {
    ## applies Deriv
    ## ignores [[
    #temp = Deriv::drule[['[[']]
    #assign('[[', list(0), envir=Deriv::drule)

    r = mapply(Deriv::Deriv, rep(list(f), length(names)), names)
    base::names(r) = names

    #assign('[[', temp, envir=Deriv::drule)
    return(r)
}

Deriv2f = function(f, names) {
    ## applies Deriv twice
    ## ignores [[
    #temp = Deriv::drule[['[[']]
    #assign('[[', list(0), envir=Deriv::drule)

    d = mapply(Deriv::Deriv, rep(list(f), length(names)), names)
    base::names(d) = names

    d2 = mapply(Deriv::Deriv, d, names)

    r = rep(list(NA), length(names))
    base::names(r) = names
    r = replicate(length(names), r, simplify=F)
    base::names(r) = names

    for (i in seq1(1, length(names))) {
        n = names[[i]]
        r[[n]][[n]] = d2[[i]]
    }

    combs = combn(names, 2)
    d2 = mapply(Deriv::Deriv, d[combs[1,]], combs[2,])

    for (i in seq1(1, ncol(combs))) {
        a = combs[1, i]
        b = combs[2, i]
        r[[a]][[b]] = d2[[i]]
        r[[b]][[a]] = d2[[i]]
    }

    #assign('[[', temp, envir=Deriv::drule)
    return(r)
}


mirrorMatrix = function(x) {
    ## transforms upper/lower diagonal matrix to full matrix
    r = x + t(x)
    diag(r) = diag(x)
    return(r)
}


## nested list helper
is_flat = function(x) !any(sapply(x, inherits, 'list'))

flatten = function(x) {
    if (!inherits(x, 'list'))
        return(list(x))
    if (is_flat(x))
        return(x)
    return(do.call(c, lapply(x, flatten) ))
}


zmin = function(x) if (length(x) == 0) 0 else min(x)
zmax = function(x) if (length(x) == 0) 0 else max(x)


lproduct = function(x) {
    ## product expands a list of lists
    if (length(x) == 0)
        return(list())
    idcs = lapply(x, seq)
    idcs = do.call(expand.grid, idcs)
    colnames(idcs) = names(x[[1]])
    r = apply(idcs, 1, function(idcs)
        mapply(function(idx, i) x[[i]][[idx]], idcs, 1:length(x), SIMPLIFY=F))
    return(r)
}


#' Integrate Alternative
#'
#' \code{integrateA} is a tolerance wrapper for \code{stats::integrate}.
#' It allows \code{integrate} to reach the maximum number of subdivisions.
#'
#' See \code{\link[stats]{integrate}}.
#'
#' @param f,lower,upper,...,subdivisions,rel.tol,abs.tol,stop.on.error,keep.xy,aux see \code{\link[stats]{integrate}}.
#'
#' @seealso \code{\link[stats]{integrate}}
#'
#' @example R/examples/integrateA.R
#'
#' @export
integrateA = function(f, lower, upper, ..., subdivisions=100L, rel.tol=.Machine$double.eps^0.25, abs.tol=rel.tol, stop.on.error=TRUE, keep.xy=FALSE, aux=NULL) {
    r = stats::integrate(f, lower, upper, ..., subdivisions=subdivisions, rel.tol=rel.tol, abs.tol=abs.tol, stop.on.error=F, keep.xy=keep.xy, aux=aux)
    if ( !(r$message %in% c('OK', 'maximum number of subdivisions reached')) ) {
        if (stop.on.error) {
            stop(r$message)
        }
    }
    return(r)
}


clusterPeak = function(x, y, maxDist) {
    ## x = row matrix of points
    ## y = corresponding vector of values (of length nrow(x))
    ##
    ## iteratively assigns points to the nearest cluster
    ## in reverse order of magnitude of y

    r = rep(0, length(y))

    dists = as.matrix(dist(x))
    idcs = seq1(1, length(y))

    rr = 1

    while (length(idcs) != 0) {
        i = which.max(y[idcs])
        idx = idcs[i] # idx of max y
        ds = dists[,idx] # distances to all other points

        cIdcs = which(ds <= maxDist) # candidates
        cIdcs = cIdcs[order(ds[cIdcs])] # sort by distance
        cSub = match(T, r[cIdcs] != 0) # first already matched
        if (is.na(cSub)) { # none matched
            r[idx] = rr
            rr = rr + 1
            next
        }

        r[idx] = r[cIdcs[cSub]] # match
        idcs = idcs[-i]
    }

    return(r)
}


#' Matrix Ordering Permutation
#'
#' \code{roworder} returns a permutation which rearranges the rows of its first argument into ascending order.
#'
#' @param x a matrix.
#' @param ... other arguments passed to \code{order}.
#'
#' @return \code{roworder} returns an integer vector.
#'
#' @seealso \code{\link[base]{order}}
#'
#' @example R/examples/roworder.R
#'
#' @export
roworder = function(x, ...) {
    cols = lapply(seq1(1, ncol(x)), function(i) x[,i])
    args = c(cols, ...)
    return(do.call(order, args))
}


#' Row Matching
#'
#' \code{rowmatch} returns a vector of the positions of (first) matches of the rows of its first argument in the rows of its second.
#'
#' \code{rowmatch} requires unique rows in \code{x} and \code{table}.
#'
#' @param x a row matrix or data.frame, the rows to be matched.
#' @param table a row matrix or data.frame, the rows to be matched against.
#' @param nomatch the value to be returned in the case when no match is found.
#' Note that it is coerced to \code{integer}.
#'
#' @return \code{rowmatch} returns an integer vector giving the position of the matching row in \code{table} for each row in \code{x}. And \code{nomatch} if there is no matching row.
#'
#' @seealso \code{\link[base]{match}}
#'
#' @example R/examples/rowmatch.R
#'
#' @export
rowmatch = function(x, table, nomatch=NA_integer_) {
    nomatch = as.integer(nomatch)

    ordx = roworder(x)
    ordt = roworder(table)
    sx = x[ordx,, drop=F]
    st = table[ordt,, drop=F]
    isx = which(duplicated(rbind(st, sx))) - nrow(st) # idcs of sorted table in sorted x
    ist = which(duplicated(rbind(sx, st))) - nrow(sx) # idcs of sorted x in sorted table

    if (any(isx < 0) || any(ist < 0))
        stop('rows in x and table shall be unique')

    ix = ordx[isx] # idcs of sorted table in x
    it = ordt[ist] # idcs of sorted x in table

    r = rep(nomatch[1], nrow(x))
    r[ix] = it
    return(r)
}

