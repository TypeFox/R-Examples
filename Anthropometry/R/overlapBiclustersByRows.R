overlapBiclustersByRows <- function(Bic,resBicluster) {

 x <- rep(0, nrow(resBicluster@RowxNumber))

 x <- x + Bic * resBicluster@RowxNumber[, Bic]

 x
}

