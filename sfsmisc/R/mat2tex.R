### Port to R and a few small improvements:
### Copyright © 2000 Martin Maechler, ETH Zurich

mat2tex <- function(x, file = "mat.tex", envir = "tabular",
                    nam.center = "l", col.center = "c",
                    append = TRUE, digits = 3, title)
{
    if(length(d.x <- dim(x)) != 2)
        stop("'x' must be a matrix like object with dim(x) of length 2")
    if(any(d.x <= 0))
        stop("'dim(x)' must be positive")
    nr.x <- d.x[1]
    nc.x <- d.x[2]
    c2ind <- (1:nc.x)[-1] # possibly empty

    ## determine if there are labels to be processed
    dn.x <- dimnames(x)
    if(has.rowlabs <- !is.null(dn.x[[1]]))        rowlabs <- dn.x[[1]]
    if(has.collabs <- !is.null(dn.x[[2]]))        collabs <- dn.x[[2]]

    ## produce column specification
    stopifnot(any(nam.center == c("l","r","c")))
    stopifnot(all(col.center %in% c("l","r","c")))
    col.center <- rep(col.center, length = nc.x)
    colspec <- "{|"
    if(has.rowlabs)
        colspec <- paste(colspec, nam.center, "||")
    colspec <- paste0(colspec, paste(col.center, "|", collapse=""), "}")
    cat(paste(sprintf("\\begin{%s}", envir), colspec, " \n"), file=file, append=append)

    span <- nc.x + if(has.rowlabs) 1 else 0
    cat(if(!missing(title)) paste("\\multicolumn{", span,
                                  "}{c}{", title, "} \\\\"),
        "\\hline \n", file = file, append = TRUE)
    ## output column labels if needed
    if(has.collabs) {
        collabline <- " "
        if(has.rowlabs)
            collabline <- paste(collabline, " \\  &")
        collabline <- paste(collabline, collabs[1])
        for(i in c2ind)
            collabline <- paste(collabline, "&", collabs[i])
        collabline <- paste(collabline, "\\\\ \\hline \\hline")
        cat(collabline, "\n", file = file, append = TRUE)
    }
    ## output matrix entries
    options(digits = digits)
    for(i in 1:nr.x) {
        thisline <-
            if(has.rowlabs)
                paste(rowlabs[i], "&", format(x[i, 1])) else format(x[i, 1])
        for(j in c2ind)
            thisline <- paste(thisline, "&", format(x[i, j]))

        thisline <- paste(thisline, "\\\\ \\hline")
        cat(paste(thisline, "\n"), file = file, append = TRUE)
    }
    cat(paste0("\\end{", envir, "}\n"), file = file, append = TRUE)
}

