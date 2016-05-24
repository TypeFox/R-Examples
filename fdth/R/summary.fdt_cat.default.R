summary.fdt_cat.default <- function (object,
                                     columns=1:6,
                                     round=2,
                                     row.names=FALSE,
                                     right=TRUE, ...)
{
  res <- cbind(object[, 1],
               round(object[, 2:6],
                     round))[columns]

  names(res) <- c('Category',
                  'f',
                  'rf',
                  'rf(%)',
                  'cf',
                  'cf(%)')[columns]

  res <- print.data.frame(res,
                          row.names=row.names,
                          right=right, ...)
}

