INV <-
function(w) {as.matrix(t(svd(w)$v %*% (t(svd(w)$u)/svd(w)$d)))}
