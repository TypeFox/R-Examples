# A "breslow" object for storing multiple survival curves on the same time scale
setClass("breslow",
  representation(
    time = "vector",
    curves = "matrix"
  )
)

setMethod("show", "breslow", function(object) {
  cat("A \"breslow\" object with", nrow(object@curves), "survival curve")
  if (nrow(object@curves)>1) cat("s")
  cat(" and", ncol(object@curves), "time points.\n")
})

setMethod("plot", "breslow", function(x, y, xlab = "time", ylab = "survival probability", ...) {
  plot(0,-1,xlim=range(x@time),ylim=0:1, ylab = ylab,xlab = xlab, ...)
  for (i in 1:nrow(x@curves)) lines(x@time, x@curves[i,], type="s", ...)
  return(invisible(NULL))
})

setMethod("as.matrix", "breslow", function(x, ...) {
  out <- x@curves
  colnames(out) <- x@time
  out
})

setGeneric("as.data.frame")
setMethod("as.data.frame", "breslow", function(x, row.names = NULL, optional = FALSE, ...) {
  subjectnames <- rownames(x@curves)
  if (is.null(subjectnames)) 
    subjectnames <- 1:nrow(x@curves)
  if (length(subjectnames) > 1) 
    col.names <- paste("survival", subjectnames, sep=".")
  else
    col.names <- "survival"
  out <- as.data.frame(t(x@curves), row.names, optional)
  colnames(out) <- make.names(col.names)
  out <- data.frame(out, time = x@time)
  out
})


setMethod("[", "breslow", function(x, i, j, ... , drop = TRUE) {
  if (missing(i) && missing(j))
    as.matrix(x)[,,drop=drop]
  else if (missing(i))
    as.matrix(x)[,j,drop=drop]
  else if (missing(j))
    as.matrix(x)[i,,drop=drop]
  else
    as.matrix(x)[i,j,drop=drop]
})

setMethod("[[", "breslow", function(x,i,j) {
  x@curves <- x@curves[i,,drop=FALSE]
  x
})

setMethod("time", "breslow", function(x, ...) {
  x@time
})


setMethod("as.list", "breslow", function(x, ...) {
  list(time = x@time, curves = x@curves)
})

setGeneric("survival", function(object, ...)  standardGeneric("survival"))
setMethod("survival", "breslow", function(object, time) {
  object[,max(which(object@time <= time))]
})
