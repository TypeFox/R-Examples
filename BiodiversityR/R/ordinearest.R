`ordinearest` <-
function (ordiplot,dist,...) {
    ord <- scores(ordiplot, display="sites",...)
    dist <- as.matrix(dist)
    diag(dist) <- Inf
    nabo <- apply(dist, 1, which.min)
    graphics::arrows(ord[,1], ord[,2], ord[nabo,1], ord[nabo,2], ...)
    invisible()
}

