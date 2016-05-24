catn <- function(...) {cat(...); cat("\n")}

pn <- function() cat("\n")

"prinV" <- function(x,after=2,before)
{ ## print vectors without []
  if (missing(before) || is.null(before))
    cat(paste(formatFix(x,after),collapse=""),"\n")
  else
    cat(paste(formatFix(x,after,before),collapse=""),"\n")
  invisible(x)
}

"prinM" <- function(x,after=2,before) {
  ## print matrices without []
  res <- x
  res <- formatFix(x, after, before)
  dim(res) <- dim(x)
  print(res,quote=FALSE)
  invisible(res)
}

"prinT" <- function(x,rownam=FALSE,colnam=FALSE) {  
  if (is.matrix(x)) {
    if (colnam) 
      cat(paste(
        if (rownam) " ", paste(dimnames(x)[[2]],collapse="\t")
      ,sep="\t"),"\n")
    for (ii in 1:dim(x)[1])
      cat(paste(
        if (rownam) dimnames(x)[[1]][ii],
	                 paste(x[ii,],collapse="\t")
      ,sep="\t"),"\n")
  }
  else {
    if (rownam) 
      cat(paste(names(x),collapse="\t"),"\n")
    cat(paste(x,collapse="\t"),"\n")
  }
  invisible(x)
}

"prinE" <- function(xsv, ...) { 
  cat(xsv, "=")
  print(as.vector(eval.parent(parse(text = xsv))), ...)
  invisible(xsv)
}

"prinL" <- function(xs, ...) {
	cat(xs, "\n")
	print(eval.parent(parse(text = xs)), ...)
  invisible(xs)
}

"prinP" <- function(xs) {
	cat(xs, "\n")
	eval.parent(parse(text = xs))
  invisible(xs)
}
NprinV <- function(x,after=2,before) {cat("\n"); prinV(x,after,before) }
NprinM <- function(x,after=2,before) {cat("\n"); prinM(x,after,before) }
NprinT <- function(x,rownam=FALSE,colnam=FALSE) {cat("\n"); prinT(x,rownam,colnam) }
NprinE <- function(xsv,...) {cat("\n"); eval.parent(substitute(prinE(xsv,...))) }
NprinL <- function(xs ,...) {cat("\n"); eval.parent(substitute(prinL(xs, ...))) }
NprinP <- function(xs)      {cat("\n"); eval.parent(substitute(prinP(xs))) }

