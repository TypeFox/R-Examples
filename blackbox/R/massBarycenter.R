massBarycenter <- function (vertices) {
    tc <- delaunayn(vertices, options = "Pp")
    pmul <- cbind(-1, diag(rep(1, ncol(vertices))))
    vb <- apply(tc, 1, function(v) {
        simplex <- vertices[v, ]
        vol <- abs(det(pmul %*% simplex))
        bary <- colMeans(simplex)
        c(vol, vol * bary)
    })
  resu <- rowSums(vb[-1, ])/sum(vb[1, ])
  return(resu)
}
