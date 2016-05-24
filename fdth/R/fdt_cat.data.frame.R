fdt_cat.data.frame <- function (x,
                                by,
                                sort=TRUE,
                                decreasing=TRUE, ...)
{
 stopifnot(is.data.frame(x))

 res <- list()

 # User do not defines a factor
 if (missing(by)) {

  x <- na.omit(x) 

  logCol <- sapply(x,
                   is.factor)

  for (i in 1:ncol(x)) {
   if (logCol[i]) {
    m <- as.data.frame(x[, i])

    fdt <- make.fdt_cat.multiple(m,
                             sort,
                             decreasing)

    tmpres <- list(table=fdt[[1]])

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
   tmpdf <- tmpdf[-pos]
   logCol <- sapply(tmpdf,
                    is.factor)

   for (j in 1:ncol(tmpdf)) {
    if (logCol[j]) {
     m <- as.data.frame(tmpdf[, j])

     fdt <- make.fdt_cat.multiple(m,
                                  sort,
                                  decreasing)
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

 class(res) <- c('fdt_cat.multiple',
                 'fdt_cat',
                 'list')

 invisible(res)
}

