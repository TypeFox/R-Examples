print.fdt.multiple <- function (x,
                                columns=1:6, 
                                round=2,
                                format.classes=FALSE,
                                pattern='%09.3e', 
                                row.names=FALSE, 
                                right=TRUE, ...)
{
  right.tmp <- as.logical(x[[1]][['breaks']]['right'])

  tnames <- names(x)

  for (i in 1:length(tnames)) {
    res <- x[tnames[i]][[tnames[i]]][['table']]

    res <- cbind(res[, 1],
                 round(res[, 2:6],
                       round))[columns]

    cat(tnames[i], '\n')

    if (format.classes) {
      tmp <- as.character(res[, 1])

      res[, 1] <- make.fdt.format.classes(tmp,
                                          right.tmp,
                                          pattern)
    }

    names(res) <- c('Class limits',
                    'f', 
                    'rf', 
                    'rf(%)', 
                    'cf',
                    'cf(%)')[columns]

    print.data.frame(res,
                     row.names=row.names,
                     right=right, ...)

    cat('\n')}
}

