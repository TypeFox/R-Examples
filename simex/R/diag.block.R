diag.block <-
function (d, n) {
if (is.list(d)) {
d.row <- sapply(d, NROW)
d.col <- sapply(d, NCOL)
d.diag <- matrix(0, nrow = sum(d.row), ncol = sum(d.col))
d.row <- c(0, cumsum(d.row))
d.col <- c(0, cumsum(d.col))
for (i in 1:length(d)) {
d.diag[(d.row[i] + 1):d.row[i + 1], (d.col[i] + 1):d.col[i + 1]] <-
as.matrix(d[[i]])
}
}
if (!is.list(d)) {
d.diag <- kronecker(diag(1, n), d)
}
return(d.diag)
}

