
isa.biclust <- function(modules) {

  if (!requireNamespace("biclust")) {
    stop("The `biclust' package is required for this")
  }

  new("Biclust", Parameters=list(seeddata=modules$seeddata,
                   rundata=modules$rundata),
      RowxNumber=modules$rows != 0,
      NumberxCol=t(modules$columns != 0), Number=ncol(modules$rows))
}
