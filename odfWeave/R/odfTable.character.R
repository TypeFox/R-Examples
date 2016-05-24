"odfTable.character" <-
function(
   x,
   horizontal = length(x) < 5,
   colnames    = names(x),  # allow alt col names (what about dimnames[[1]]?)
   name = paste("Table", floor(runif(1) * 1000), sep = ""),
   styles = NULL,
   ...)
{
   if(!is.null(colnames)) colnames <- odfTranslate(colnames, toR = FALSE)
   xMat <- if(horizontal) as.matrix(t(x)) else as.matrix(x)
   colTypes <- apply(xMat, 2, odfDataType)

   theDots <- list(...)
   if(!any(names(theDots) == "justify")) theDots$justify <- "none"
   if(!any(names(theDots) == "trim")) theDots$trim <- TRUE
   
   args <- list(x = xMat)
   args <- c(args, theDots)
   xChar <- do.call("format", args)

   if(is.null(styles)) styles <- tableStyles(xChar, useRowNames = FALSE, colnames)

   tbleText <- odfTableGen(xChar, colTypes, header = colnames, tableName = name, styles)
   structure(tbleText, class = "odfTable")
}


"odfTable.factor" <-
function(
   x,
   horizontal = length(x) < 5,
   colnames    = names(x),  # allow alt col names (what about dimnames[[1]]?)
   name = paste("Table", floor(runif(1) * 1000), sep = ""),
   styles = NULL,
   ...)
{
   if(!is.null(colnames)) colnames <- odfTranslate(colnames, toR = FALSE)
   odfTable.character(as.character(x), horizontal = horizontal, colnames = colnames, name = name, styles = styles, ...)
}

