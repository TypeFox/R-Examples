print.coef <- function(x, dec=3){
    # suspend warnings for this function:
    w <- options()$warn; options(warn=-1)
    # check coefficient names in the list
    if (!all(nchar(names(x))) || is.null(names(x))){  
        # note: all(integer(0))=T (i.e. no-names at all would still result 
        # F to the first part of the if statement, hence the second part
        ind <- if (is.null(names(x))) seq(x) else !nchar(names(x))
        names(x)[ind] <- paste('X', seq(x)[ind], sep='')
    }  
    k <- unlist(lapply(x, length))
    ind <- order(k, decreasing=T)
    x <- x[ind]; k <- k[ind]
    # set up output table:
    tab <- round(data.frame(x[c(1,1)]), dec)
    names(tab)[2] <- paste(names(x)[1], 'c', sep='.')
    tab[1] <- rownames(tab); rownames(tab) <- NULL
    for(i in seq(x[-1])){ # run for the remaining coefficients in the list
        tab[ncol(tab)+1] <- NA 
        if(!all(names(x[[i]])==names(x[[i+1]]))) {
            names(tab)[ncol(tab)] <- ' ' 
            # this column name will be adjusted at repeated passes, so 
            # it needs one more correction at the end (see below)
            tab[[names(x)[i+1]]] <- NA
            tab[[ncol(tab)]][seq(k[i+1])] <- names(x[[i+1]])
            tab[ncol(tab)+1] <- NA            
        }  
        tab[[ncol(tab)]][seq(k[i+1])] <- round(x[[i+1]], dec)
        names(tab)[ncol(tab)] <- paste(names(x)[i+1], 'c', sep='.')
    }
    tab[is.na(tab)] <- '  '
    # final correction of empty column names:
    names(tab)[substring(names(tab), 1, 1)==' '] <- '  '
    options(warn=w)
    print(tab)  
}
