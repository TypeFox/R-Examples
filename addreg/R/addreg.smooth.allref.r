addreg.smooth.allref <- function(object, data = environment(object), mono, family,
                                    addreg.smooth.spec, num.knots) {
    t <- if(missing(data))
        terms(object)
    else terms(object, data = data)

    termlist <- attr(t, "term.labels")
    nvar <- length(termlist)
    smoothlist <- sapply(addreg.smooth.spec$smooth.spec,"[[","term")
    smoothtype <- sapply(addreg.smooth.spec$smooth.spec,class)
    names(smoothtype) <- smoothlist
    nsmvar <- length(smoothlist)
    
    if(length(num.knots) != nsmvar)
        stop(gettextf("num.knots has length %d should equal %d (number of smooth terms)",
                length(num.knots), nsmvar), domain = NA)
    num.knots <- as.vector(num.knots, mode = "integer")
    names(num.knots) <- smoothlist
    
    if(missing(mono)) mono <- rep(FALSE, nvar)
    if(is.null(mono)) mono <- rep(FALSE, nvar)
    monotonic <- rep(FALSE, nvar)
    names(monotonic) <- termlist
    monotonic[mono] <- TRUE
	names(monotonic) <- termlist
    
    allref <- list()
    for(smth in smoothlist) {
        allref[[smth]] <- list()
        if(smoothtype[smth] == "Iso.smooth") {
            allref[[smth]][[1]] <- 1
        } else if(smoothtype[smth] == "B.smooth") {
            if(monotonic[smth]) allref[[smth]][[1]] <- 1:(num.knots[smth]+3)
            else {
                if(family$family != "binomial") 
                    allref[[smth]] <- as.list(1:(num.knots[smth]+3))
                else
                    allref[[smth]] <- combinat::permn(1:(num.knots[smth]+3))
            }
        } else
            stop("smooth type not recognized. Only B() and Iso() are supported by addreg.smooth")
    }
    
    list(allref = allref, terms = t, data = data, monotonic = monotonic)
}