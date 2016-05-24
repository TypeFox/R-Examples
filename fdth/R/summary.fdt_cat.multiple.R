summary.fdt_cat.multiple <- function (object,
                                      columns=1:6,
                                      round=2,
                                      row.names=FALSE,
                                      right=TRUE, ...)
{
  tnames <- names(object)

  for (i in 1:length(tnames)) {
    res <- object[[tnames[i]]]

    res <- cbind(res[[1]][, 1],
                 round(res[[1]][, 2:6],
                       round))[columns]

    cat(tnames[i], '\n')

    names(res) <- c('Category',
                    'f', 
                    'rf', 
                    'rf(%)',
                    'cf',
                    'cf(%)')[columns]

    print.data.frame(res,
                     row.names=row.names, 
                     right=right, ...)

    cat('\n')
  }
}

