
isep <- function(x)
    paste(paste(x[-length(x)], "&", collapse = " "), 
          x[length(x)], collapse = " ")

caption <- function(xname, label, caption, pkg = NULL) {
    RET <- paste("\\caption{\\Robject{", xname, "} data", 
                 sep = "", collapse = "")
    if (!is.null(pkg))
        RET <- paste(RET, " (package \\Rpackage{", pkg, "})", 
                     sep = "", collapse = "")
    RET <- paste(RET, ". ", caption, sep = "", collapse = "")
    RET <- paste(RET, paste("\\label{", label, "}}", 
                 sep = "", collapse = ""))
    return(RET)
}


HSAURtable <- function(object, ...)
    UseMethod("HSAURtable")

HSAURtable.data.frame <- function(object, xname = deparse(substitute(object)),
                                  pkg = NULL, nrows = NULL,...) {

    digits <- 0:6
    table <- matrix("0", nrow = nrow(object), ncol = ncol(object))
    xcc <- object[complete.cases(object),]
    for (i in 1:ncol(object)) {
        if (is.numeric(xcc[[i]])) {
            d <- min(which(sapply(digits, 
                function(d) 
                    max(abs(xcc[[i]] - round(xcc[[i]], d))) < 
                            sqrt(.Machine$double.eps))))
            table[,i] <- formatC(object[[i]], digits = digits[d], format = "f")
        } else {
            table[,i] <- as.character(object[[i]])
        }
    }
    if (!is.null(nrows)) table <- rbind(table[1:nrows,,drop = FALSE], "$\\vdots$")

    RET <- list(xname = xname,
                pkg = pkg, 
                varnames = colnames(object),
                rownames = rownames(object),
                data = table)
    class(RET) <- "dftab"
    return(RET)
}

HSAURtable.table <- function(object, xname = deparse(substitute(object)),
                             pkg = NULL,...) {

    xtab <- matrix(as.character(object), nrow = nrow(object), 
                   ncol = ncol(object))
    RET <- list(xname = xname,
                pkg = pkg,
                varnames = names(dimnames(object)),
                data = rbind(c(" ", dimnames(object)[[2]]),
                             cbind(dimnames(object)[[1]], xtab)))
    class(RET) <- "tabtab"
    return(RET)
}

toLatex.tabtab <- function(object, caption = NULL, label = NULL, 
                           topcaption = TRUE, index = TRUE, ...) {

    RET <- c()
    nc <- ncol(object$data)

    if (index)
        RET <- c(RET, paste("\\index{", object$xname, " data@\\Robject{",
                            object$xname, "} data}", sep = ""))

    RET <- c(RET, "\\begin{center}")

    RET <- c(RET, paste("\\begin{longtable}",
        paste("{", paste(rep("r", nc + 1), collapse = ""), "}")))
    if (topcaption)
        RET <- c(RET, caption(object$xname, label, caption, object$pkg),
                 "\\\\")
    RET <- c(RET, paste(" & & \\multicolumn{", nc - 1, "}{c}{\\Robject{", 
              object$varnames[2], "}} \\\\", collapse = ""))
    object$data <- cbind(c(paste("\\Robject{", object$varnames[1], "}", 
                                 collapse = ""), 
                           rep(" ", nrow(object$data) - 1)), object$data)
    RET <- c(RET,  apply(object$data, 1, function(x) paste(isep(x), "\\\\"))) 
    if (!topcaption)
        RET <- c(RET, caption(object$xname, label, caption, object$pkg))
    RET <- c(RET, "\\end{longtable}")
    RET <- c(RET, "\\end{center}")
    class(RET) <- "Latex"
    return(RET)
}


toLatex.dftab <- function(object, pcol = 1, caption = NULL, 
    label = NULL, rownames = FALSE, topcaption = TRUE, index = TRUE, ...) {
    
    nc <- ncol(object$data)

    if (pcol > 1) {
        nr <- ceiling(nrow(object$data) / pcol)
        object$data <- rbind(object$data, matrix(" ", 
            nrow = nr * pcol - nrow(object$data), 
            ncol = nc))
        d <- NULL
        for (i in 1:pcol)
            d <- cbind(d, object$data[((i - 1) * nr + 1):(i * nr),])
        object$data <- d       
    }

    RET <- c()

    if (index)
        RET <- c(RET, paste("\\index{", object$xname, " data@\\Robject{",
                            object$xname, "} data}", sep = ""))

    RET <- c(RET, "\\begin{center}")
    if (rownames) 
        RET <- c(RET, 
            paste("\\begin{longtable}{l", paste(rep(paste(rep("r", nc), 
                                                          collapse = ""), pcol), 
                                                collapse = "|"), "}", 
                  collapse = ""))
    else 
        RET <- c(RET, 
            paste("\\begin{longtable}{", paste(rep(paste(rep("r", nc), 
                                                         collapse = ""), pcol), 
                                              collapse = "|"), "}", 
                  collapse = ""))
    if (topcaption)
        RET <- c(RET, caption(object$xname, label, caption, object$pkg),
                 "\\\\")
    RET <- c(RET, "\\hline")
    vn <- rep(object$varnames, pcol)
    vn <- paste(paste("\\Robject{", vn, sep = ""), "}", sep = "")
    if (rownames) {
        RET <- c(RET, paste("  & ", isep(vn), "\\\\ \\hline"))
        RET <- c(RET, "\\endfirsthead")
        RET <- c(RET, paste("\\caption[]{\\Robject{", object$xname, 
                            "} data (continued).} \\\\", 
                 sep = "", collapse = ""))
        RET <- c(RET, "\\hline")
        RET <- c(RET, paste("  & ", isep(vn), "\\\\ \\hline"))
        RET <- c(RET, "\\endhead")
        for (i in 1:nrow(object$data))
            RET <- c(RET, paste(object$rownames[i], "  & ", 
                                isep(object$data[i,]), "\\\\"))
    } else {
        RET <- c(RET, paste(isep(vn), "\\\\ \\hline"))
        RET <- c(RET, "\\endfirsthead")
        RET <- c(RET, paste("\\caption[]{\\Robject{", object$xname, 
                            "} data (continued).} \\\\", sep = "", collapse = ""))
        RET <- c(RET, "\\hline")
        RET <- c(RET, paste(isep(vn), "\\\\ \\hline"))
        RET <- c(RET, "\\endhead")
        RET <- c(RET, 
            apply(object$data, 1, function(x) paste(isep(x), "\\\\")))
    }
    RET <- c(RET, "\\hline")
    if (!topcaption)
        RET <- c(RET, caption(object$xname, label, caption, object$pkg))
    RET <- c(RET, "\\end{longtable}")
    RET <- c(RET, "\\end{center}")
    class(RET) <- "Latex"
    return(RET)
}

