matches <-
function (subs, mat) {
    m1 <- subs
    m2 <- mat[, -ncol(mat)]
    stopifnot(ncol(m1) == ncol(m2), typeof(m1) == typeof(m2)) 
    m1 <- apply(m1, 1, paste, collapse = "")
    m2 <- apply(m2, 1, paste, collapse = "")
    mat[match(m1, m2), ]
}
