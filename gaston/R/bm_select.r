select.snps <- function(x, condition) {
  if(!is(x, "bed.matrix")) stop("x is not a bed.matrix")
  w <- eval(substitute(condition), x@snps, parent.frame())
  x[,w]
}

select.inds <- function(x, condition) {
  if(!is(x, "bed.matrix")) stop("x is not a bed.matrix")
  w <- eval(substitute(condition), x@ped, parent.frame())
  x[w,]
}

