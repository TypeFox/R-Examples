print.fdt_cat.default <- function (x,
                                   columns=1:6,
                                   round=2,
                                   row.names=FALSE,
                                   right=TRUE, ...)
{
  res <- cbind(x[, 1],
               round(x[, 2:6],
                     round))[columns]

  names(res) <- c('Category',
                  'f',
                  'rf',
                  'rf(%)',
                  'cf',
                  'cf(%)')[columns]

  print.data.frame(res, 
                   row.names=row.names, 
                   right=right, ...)
}

