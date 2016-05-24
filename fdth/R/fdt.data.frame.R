fdt.data.frame <- function (x,
                            k,
                            by,
                            breaks=c('Sturges', 'Scott', 'FD'),
                            right=FALSE, ...)
{
  stopifnot(is.data.frame(x))

  res <- list()

  # User do not defines a factor
  if (missing(by)) {
    logCol <- sapply(x,
                     is.numeric)

    for (i in 1:ncol(x)) {
      if (logCol[i]) {
        m <- as.matrix(x[, i])

        fdt <- make.fdt.multiple(m,
                                 k,
                                 breaks,
                                 right)

        tmpres <- list(table=fdt[[1]],
                       breaks=fdt[[2]])

        res <- c(res,
                 list(tmpres))
      }
    }

    valCol     <- logCol[logCol]

    names(res) <- names(valCol)
  }

  # User defines one factor
  else {
    nameF   <- character()
    nameY   <- character()
    namesdf <- names(x)
    pos     <- which(namesdf == by)

    stopifnot(is.factor((x[[pos]])))

    numF <- table(x[[pos]])
    for (i in 1:length(numF)) {
      tmpdf  <- subset(x,
                       x[[pos]] == names(numF[i]))

      logCol <- sapply(tmpdf,
                       is.numeric)

      for (j in 1:ncol(tmpdf)) {
        if (logCol[j]) {
          m <- as.matrix(tmpdf[, j])

          fdt <- make.fdt.multiple(m,
                                   k,
                                   breaks,
                                   right)
          newFY  <- list(fdt)
          nameF  <- names(numF[i])
          nameY  <- names(logCol[j])
          nameFY <- paste(nameF,
                          '.',
                          nameY,
                          sep="")

          names(newFY) <- sub(' +$',
                              '',
                              nameFY)
          res <- c(res,
                   newFY)
        }
      }
    }
  }

  class(res) <- c('fdt.multiple',
                  'fdt',
                  'list')
  invisible(res)
} 
