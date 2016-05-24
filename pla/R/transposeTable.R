.transposeTable <- function(TableChar, firstName, RowNames, ColumnNames, 
                            fun, log) {
    Nsamples  <- length(unique(ColumnNames))
    Nsteps    <- dim(TableChar)[2]
    if (Nsamples > 1) {
        if (firstName == "Dose") {
            Nsteps <- dim(TableChar)[1]/Nsamples
            Table <- matrix(as.numeric(t(TableChar)),
                            nrow = dim(TableChar)[2])
            Table <- fun(Table)
            if (log == TRUE) 
                Table <- log(Table)
            else if (log != FALSE) 
                Table <- log(Table)/log(log)
            dimnames(Table)[[2]] <- paste(rep(unique(ColumnNames),
                                              rep(Nsteps,
                                                  Nsamples)),
                                          rep(1:Nsteps, Nsamples),
                                          RowNames,
                                          sep = ":")
        }
        else {
            Table <- matrix(as.numeric(as.matrix(TableChar)),
                            ncol = Nsteps * Nsamples)
            Table <- fun(Table)
            if (log == TRUE) 
                Table <- log(Table)
            else if (log != FALSE) 
                Table <- log(Table)/log(log)
            nmsDoses <- dimnames(TableChar)[[2]]
            nmsDoses <- unlist(lapply(strsplit(nmsDoses, ":"),
                                      FUN = function (x)
                                      ifelse(length(x) > 1, x[2],
                                             substr(x, 2, nchar(x))
                                             )))
            dimnames(Table)[[2]] <- paste(rep(unique(ColumnNames),
                                              Nsteps),
                                          rep(1:Nsteps,
                                              rep(Nsamples,
                                                  Nsteps)),
                                          rep(nmsDoses,
                                              rep(Nsamples,
                                                  Nsteps)),
                                          sep = ":")
        }
    }
    else {
        Table <- matrix(as.numeric(as.matrix(TableChar)),
                        ncol = dim(TableChar)[2])
        Table <- fun(Table)
        if (log) 
            Table <- log(Table)
        if (all(dim(Table) == dim(TableChar))) {
            dimnames(Table) <- dimnames(TableChar)
            nms <- dimnames(Table)[[2]]
            nms <- lapply(strsplit(nms, "\\."),
                          FUN = function(x) paste(x,
                              collapse = ":"))
            dimnames(Table)[[2]] <- nms
        }
    }
    return(Table)
}
