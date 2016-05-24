
forODFTable <- function(x, digits = 1, ...)
{
    if(class(x) != "CrossTable"){
        msg <- sprintf(gettext("'%s' should be of class 'CrossTable'.",
                               domain = "R-descr"), deparse(substitute(x)))
        stop(msg)
    }

    if(x$format == "SPSS")
        hdd <- 100
    else
        hdd <- 1
    nr <- dim(x$tab)[1]
    nc <- dim(x$tab)[2]

    tab <- format(x$tab, ...)
    if(x$expected == TRUE){
        xex <- outer(x$rs, x$cs, "*")
        xex <- xex / x$gt
        xx <- format(round(xex, digits), trim = TRUE, ...)
        tab <- paste(tab, xx, sep = "<text:line-break/>")
        tab <- matrix(tab, nrow = length(x$rs), ncol = length(x$cs))
    }
    if(x$prop.chisq){
        xx <- ((x$CST$expected - x$tab) ^ 2) / x$CST$expected
        xx <- format(round(xx, digits), trim = TRUE, ...)
        tab <- paste(tab, xx, sep = "<text:line-break/>")
        tab <- matrix(tab, nrow = length(x$rs), ncol = length(x$cs))
    }
    if(!is.na(x$prop.row[1])){
        xx <- format(round(x$prop.row * hdd, digits), trim = TRUE, ...)
        if(hdd == 100)
            xx <- matrix(paste(xx, "%", sep = ""), nrow = nr, ncol = nc)
        tab <- paste(tab, xx, sep = "<text:line-break/>")
        tab <- matrix(tab, nrow = length(x$rs), ncol = length(x$cs))
    }
    if(!is.na(x$prop.col[1])){
        xx <- format(round(x$prop.col * hdd, digits), trim = TRUE, ...)
        if(hdd == 100)
            xx <- matrix(paste(xx, "%", sep = ""), nrow = nr, ncol = nc)
        tab <- paste(tab, xx, sep = "<text:line-break/>")
        tab <- matrix(tab, nrow = length(x$rs), ncol = length(x$cs))
    }
    if(!is.na(x$prop.tbl[1])){
        xx <- format(round(x$prop.tbl * hdd, digits), trim = TRUE, ...)
        if(hdd == 100)
            xx <- matrix(paste(xx, "%", sep = ""), nrow = nr, ncol = nc)
        tab <- paste(tab, xx, sep = "<text:line-break/>")
        tab <- matrix(tab, nrow = length(x$rs), ncol = length(x$cs))
    }
    if(!is.na(x$resid) && x$resid == TRUE && x$expected == TRUE){
        xx <- x$tab - xex
        xx <- format(round(xx, digits), trim = TRUE, ...)
        tab <- paste(tab, xx, sep = "<text:line-break/>")
        tab <- matrix(tab, nrow = length(x$rs), ncol = length(x$cs))
    }
    if(!is.na(x$sresid) && x$sresid == TRUE && x$expected == TRUE){
        xx <- x$CST$residual
        xx <- format(round(xx, digits), trim = TRUE, ...)
        tab <- paste(tab, xx, sep = "<text:line-break/>")
        tab <- matrix(tab, nrow = length(x$rs), ncol = length(x$cs))
    }
    if(!is.na(x$asr[1])){
        xx <- format(round(x$asr, digits), trim = TRUE, ...)
        tab <- paste(tab, xx, sep = "<text:line-break/>")
        tab <- matrix(tab, nrow = length(x$rs), ncol = length(x$cs))
    }
    tab <- cbind(tab, x$rs)
    tab <- rbind(tab, c(x$cs, x$gt))
    rownames(tab)[dim(tab)[1]] <- gettext("Total", domain = "R-descr")
    colnames(tab)[dim(tab)[2]] <- gettext("Total", domain = "R-descr")

    tab
}

