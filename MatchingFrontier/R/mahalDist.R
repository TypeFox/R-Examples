mahalDist <-
function(Tnms, Cnms, inv.cov, dataset) {
    stopifnot(!is.null(dimnames(inv.cov)[[1]]), dim(inv.cov)[1] >
              1, all.equal(dimnames(inv.cov)[[1]], dimnames(inv.cov)[[2]]),
              all(dimnames(inv.cov)[[1]] %in% names(dataset)))
    covars <- dimnames(inv.cov)[[1]]
    xdiffs <- as.matrix(dataset[Tnms, covars])
    xdiffs <- xdiffs - as.matrix(dataset[Cnms, covars])
    rowSums((xdiffs %*% inv.cov) * xdiffs)
}
