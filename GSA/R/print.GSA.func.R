print.GSA.func <- function(x, ...) {
  cat("Call:\n")
  dput(x$call)
o1=rank(x$scores)>(length(x$scores)-10)
oo1=order(-x$scores[o1])
geneset.names=x$geneset.names
if(is.null(geneset.names)){geneset.names=as.character(1:length(x$scores))}
  mat1=cbind((1:length(x$scores))[o1],x$geneset.names[o1], round(x$scores[o1],4))[oo1,]

o2=rank(x$scores)<10
oo2=order(x$scores[o2])
geneset.names=x$geneset.names
if(is.null(geneset.names)){geneset.names=1:length(x$scores)}
  mat2=cbind((1:length(x$scores))[o2],x$geneset.names[o2], round(x$scores[o2],4))[oo2,]


  print(mat1, quote = FALSE)
print("")
print("")

  print(mat2, quote = FALSE)

  invisible()
}

