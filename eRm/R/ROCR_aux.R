## ---------------------------------------------------------------------------
## Dealing with argument lists, especially '...'
## ---------------------------------------------------------------------------

## return list of selected arguments, skipping those that
## are not present in arglist
.select.args <- function( arglist, args.to.select, complement=FALSE) {
    match.bool <- names(arglist) %in% args.to.select
    if (complement==TRUE) match.bool <- !match.bool
    return( arglist[ match.bool] )
}

## return arguments in arglist which match prefix, with prefix removed
## ASSUMPTION: prefix is separated from rest by a '.'; this is removed along
## with the prefix
.select.prefix <- function( arglist, prefixes, complement=FALSE ) {
    match.expr <- paste(paste('(^',prefixes,'\\.)',sep=""),collapse='|')
    match.bool <- (1:length(arglist)) %in% grep( match.expr, names(arglist) )
    if (complement==TRUE) match.bool <- !match.bool
    arglist <- arglist[ match.bool]
    names(arglist) <- sub( match.expr, '', names(arglist))
    
    return( arglist )
}

.garg <- function( arglist, arg, i=1) {
    if (is.list(arglist[[arg]])) arglist[[ arg ]][[i]]
    else arglist[[ arg ]]
}

.sarg <- function( arglist, ...) {
    ll <- list(...)
    for (argname in names(ll) ) {
        arglist[[ argname ]] <- ll[[ argname ]]
    }
    return(arglist)
}

.farg <- function( arglist, ...) {
    ll <- list(...)
    for (argname in names(ll) ) {
        if (length(arglist[[argname]])==0)
          arglist[[ argname ]] <- ll[[ argname ]]
    }
    return(arglist)
}

.slice.run <- function( arglist, runi=1) {
    r <- lapply( names(arglist), function(name) .garg( arglist, name, runi))
    names(r) <- names(arglist)
    r
}

## ---------------------------------------------------------------------------
## Line segments
## ---------------------------------------------------------------------------

.construct.linefunct <- function( x1, y1, x2, y2) {
    if (x1==x2) {
        stop("Cannot construct a function from data.")
    }

    lf <- eval(parse(text=paste("function(x) {",
        "m <- (",y2,"-",y1,") / (",x2,"-",x1,");",
        "c <- ",y1," - m * ",x1,";",
        "return( m * x + c)}",sep=" ")))
    lf
}

.intersection.point <- function( f, g ) {
    ## if lines are parallel, no intersection point
    if (f(1)-f(0) == g(1)-g(0)) {
        return( c(Inf,Inf) )
    }

    ## otherwise, choose search interval
    imin <- -1
    imax <- 1
    while (sign(f(imin)-g(imin)) == sign(f(imax)-g(imax))) {
        imin <- 2*imin
        imax <- 2*imax
    }
    h <- function(x) { f(x) - g(x) }

    intersect.x <- uniroot( h, interval=c(imin-1,imax+1) )$root
    intersect.y <- f( intersect.x )
    return( c(intersect.x, intersect.y ))
}
