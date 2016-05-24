## parms2plot <- function(parnames, parms, regex, random){
##     ## Old function
##     re.leaf <- "[]\\[_(),[:digit:][:space:]]"
##     stems <- unique(gsub(re.leaf, "", parnames)) # unique parameter "groups"
##     if (is.null(parms) && is.null(regex)){
##         ## Create a regular expression that selects all in 'parnames'.
##         re <- paste("^", stems, re.leaf, "*$", sep="")
##     } else {
##         if (!is.null(parms)){
##             ## 'parms' will eventually be used as a regular expression,
##             ## so all special characters currently in 'parms' need to
##             ## have backslashes added to them
##             parms <- gsub("(\\[|\\]|\\.|\\+|\\*|\\(|\\))", "\\\\\\1", parms)
##             parms <- paste("^", parms, re.leaf, "*$", sep="")
##         }
##         re <- c(parms, regex)
##     }
##     parlist <- lapply(re, function(r, p) p[grep(r, p)], p=parnames)
##     if (!is.null(random)){
##         random <- rep(random, length=length(re))
##         random[is.na(random)] <- length(parnames)
##         for (i in seq(along=parlist)){
##             x <- parlist[[i]]
##             r <- random[i]
##             if (length(x)>r)
##                 parlist[[i]] <- x[sort(sample(seq(along=x), r))]
##         }
##     }
##     parnames <- unlist(parlist)
##     return(parnames)
## }

parms2plot <- function(parnames, parms, regex, random, leaf.marker="[\\[_]", do.unlist=TRUE){
    addBackslash <- function(x){
        ## helper function
        ## adds a backslash to special characters so a string
        ## can be converted to a regex
        ## if (is.null(x)) return(NULL)
        gsub("(\\[|\\]|\\.|\\+|\\*|\\(|\\))", "\\\\\\1", x)
    }
    ## Replace numbers in parms with corresponding parameter names
    if (!is.null(parms)){
        param.num <- suppressWarnings(as.numeric(parms))
        parnumidx <- which(!is.na(param.num))
        parms[parnumidx] <- parnames[param.num[parnumidx]]
    }
    re.leaf <- paste(leaf.marker, ".*$", sep="")
    plot.all <- is.null(parms) && is.null(regex)
    re <- parlist.names <- NULL
    if (plot.all){
        ## Can't just return(parnames) if all parameters are desired
        ## because random may be specified.
        has.leaf <- grepl(re.leaf, parnames)
        if (any(!has.leaf)){
            parlist.names <- c(parlist.names, parnames[!has.leaf])
            re <- c(re, paste("^", addBackslash(parnames[!has.leaf]), "$", sep="")) # exact matches
        }
        if (any(has.leaf)){
            parlist.names <- c(parlist.names, unique(gsub(re.leaf, "", parnames[has.leaf])))
            re <- c(re, paste("^", addBackslash(unique(gsub(re.leaf, "", parnames[has.leaf]))), re.leaf, sep=""))
        }
    } else {
        parlist.names <- c(parms, regex)
        if (!is.null(parms)){
            exact.idx <- parms %in% parnames
            exact.matches <- parms[exact.idx]
            re <- parms
            re[exact.idx] <- paste("^", addBackslash(exact.matches), "$", sep="")
            re[!exact.idx] <- paste("^", addBackslash(parms[!exact.idx]), re.leaf, sep="")
        }
        re <- c(re, regex)
    }
    parlist <- lapply(re, function(r, p) grep(r, p, value=TRUE), p=parnames)
    if (!is.null(random)){
        random <- rep(random, length=length(re))
        random[is.na(random)] <- length(parnames) # select all parameters for groups with random=NA
        for (i in seq(along=parlist)){
            x <- parlist[[i]]
            r <- random[i]
            if (length(x)>r)
                parlist[[i]] <- x[sort(sample(seq(along=x), r))]
        }
    }
    if (do.unlist) parlist <- unlist(parlist)
    else names(parlist) <- parlist.names
    return(parlist)
}
