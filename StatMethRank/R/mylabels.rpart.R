mylabels.rpart <- function (object, digits = 4, minlength = 1L, pretty, collapse = TRUE) 
{
    if (missing(minlength) && !missing(pretty)) {
        minlength <- if (is.null(pretty)) 
            1L
        else if (is.logical(pretty)) {
            if (pretty) 
                4L
            else 0L
        }
        else 0L
    }
    ff <- object$frame
    n <- nrow(ff)
    if (n == 1L) 
        return("root")
    is.leaf <- (ff$var == "<leaf>")
    whichrow <- !is.leaf 
    vnames <- ff$var[whichrow]
    
    index <- seq(length(vnames))
    irow <- index[c(whichrow, FALSE)]
    irow <- seq(length(vnames))
    ncat <- object$splits[irow, 2L]
    lsplit <- rsplit <- character(length(irow))
    if (any(ncat < 2L)) {
        jrow <- irow[ncat < 2L]
        cutpoint <- myformatg(object$splits[jrow, 3L], digits)
        temp1 <- (ifelse(ncat < 0, "< ", ">="))[ncat < 2L]
        temp2 <- (ifelse(ncat < 0, ">=", "< "))[ncat < 2L]
        lsplit[ncat < 2L] <- paste0(temp1, cutpoint)
        rsplit[ncat < 2L] <- paste0(temp2, cutpoint)
    }
    if (any(ncat > 1L)) {
        xlevels <- attr(object, "xlevels")
        jrow <- seq_along(ncat)[ncat > 1L]
        crow <- object$splits[irow[ncat > 1L], 3L]
        cindex <- (match(vnames, names(xlevels)))[ncat > 1L]
        if (minlength == 1L) {
            if (any(ncat > 52L))
                warning("more than 52 levels in a predicting factor, truncated for printout", 
                  domain = NA)
            xlevels <- lapply(xlevels, function(z) c(letters, 
                LETTERS)[pmin(seq_along(z), 52L)])
        }
        else if (minlength > 1L) 
            xlevels <- lapply(xlevels, abbreviate, minlength)
        for (i in seq_along(jrow)) {
            j <- jrow[i]
            splits <- object$csplit[crow[i], ]
            cl <- if (minlength == 1L) 
                ""
            else ","
            lsplit[j] <- paste((xlevels[[cindex[i]]])[splits == 
                1L], collapse = cl)
            rsplit[j] <- paste((xlevels[[cindex[i]]])[splits == 
                3L], collapse = cl)
        }
    }
    if (!collapse) {
        ltemp <- rtemp <- rep("<leaf>", n)
        ltemp[whichrow] <- lsplit
        rtemp[whichrow] <- rsplit
        return(cbind(ltemp, rtemp))
    }
    lsplit <- paste0(ifelse(ncat < 2L, "", "="), lsplit)
    rsplit <- paste0(ifelse(ncat < 2L, "", "="), rsplit)
    varname <- (as.character(vnames))
    node <- as.numeric(row.names(ff))
    parent <- match(node%/%2L, node[whichrow])
    odd <- (as.logical(node%%2L))
    labels <- character(n)
    labels[odd] <- paste0(varname[parent[odd]], rsplit[parent[odd]])
    labels[!odd] <- paste0(varname[parent[!odd]], lsplit[parent[!odd]])
    labels[1L] <- "root"
    labels
}