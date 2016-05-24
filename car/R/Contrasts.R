# last modified 2 Dec 2002 by J. Fox
# all of these functions are adapted from functions in the R base package

contr.Treatment <- function (n, base = 1, contrasts = TRUE) {
    if (is.numeric(n) && length(n) == 1) 
        levs <- 1:n
    else {
        levs <- n
        n <- length(n)
    }
    lev.opt <- getOption("decorate.contrasts")
    pre <- if (is.null(lev.opt)) "[" else lev.opt[1]
    suf <- if (is.null(lev.opt)) "]" else lev.opt[2]
    dec <- getOption("decorate.contr.Treatment")
    dec <- if (!contrasts) ""
           else if (is.null(dec)) "T." 
           else dec
    contr.names <- paste(pre, dec, levs, suf, sep="")
    contr <- array(0, c(n, n), list(levs, contr.names))
    diag(contr) <- 1
    if (contrasts) {
        if (n < 2) 
            stop(paste("Contrasts not defined for", n - 1, "degrees of freedom"))
        if (base < 1 | base > n) 
            stop("Baseline group number out of range")
        contr <- contr[, -base, drop = FALSE]
    }
    contr
}

contr.Sum <- function (n, contrasts = TRUE) 
{
    if (length(n) <= 1) {
        if (is.numeric(n) && length(n) == 1 && n > 1) 
            levels <- 1:n
        else stop("Not enough degrees of freedom to define contrasts")
    }
    else levels <- n
    lenglev <- length(levels)
    lev.opt <- getOption("decorate.contrasts")
    pre <- if (is.null(lev.opt)) "[" else lev.opt[1]
    suf <- if (is.null(lev.opt)) "]" else lev.opt[2]
    dec <- getOption("decorate.contr.Sum")
    dec <- if (!contrasts) ""
           else if (is.null(dec)) "S." 
           else dec
    show.lev <- getOption("contr.Sum.show.levels")
    contr.names <- if ((is.null(show.lev)) || show.lev) paste(pre, dec, levels, suf, sep="")
    if (contrasts) {
        cont <- array(0, c(lenglev, lenglev - 1), list(levels, 
            contr.names[-lenglev]))
        cont[col(cont) == row(cont)] <- 1
        cont[lenglev, ] <- -1
    }
    else {
        cont <- array(0, c(lenglev, lenglev), list(levels,
            contr.names))
        cont[col(cont) == row(cont)] <- 1
    }
    cont
}


contr.Helmert <- function (n, contrasts = TRUE) 
{
    if (length(n) <= 1) {
        if (is.numeric(n) && length(n) == 1 && n > 1) 
            levels <- 1:n
        else stop("contrasts are not defined for 0 degrees of freedom")
    }
    else levels <- n
    lenglev <- length(levels)
    lev.opt <- getOption("decorate.contrasts")
    pre <- if (is.null(lev.opt)) "[" else lev.opt[1]
    suf <- if (is.null(lev.opt)) "]" else lev.opt[2]
    dec <- getOption("decorate.contr.Helmert")
    dec <- if (!contrasts) ""
           else if (is.null(dec)) "H." 
           else dec
    nms <- if (contrasts) 1:lenglev else levels
    contr.names <- paste(pre, dec, nms, suf, sep="")
    if (contrasts) {
        cont <- array(-1, c(lenglev, lenglev - 1), list(levels, 
            contr.names[-lenglev]))
        cont[col(cont) <= row(cont) - 2] <- 0
        cont[col(cont) == row(cont) - 1] <- 1:(lenglev - 1)
    }
    else {
        cont <- array(0, c(lenglev, lenglev), list(levels, contr.names))
        cont[col(cont) == row(cont)] <- 1
    }
    cont
}
