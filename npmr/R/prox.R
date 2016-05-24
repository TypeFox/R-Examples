prox <-
function(B, threshold, group) {
    B. = B
    for (g in 1:max(group)) {
        SVD = svd(B[group == g, ])
        D = (SVD$d - threshold)*(SVD$d - threshold > 0)
        B.[group == g, ] = SVD$u %*% (D * t(SVD$v))
    }
    B.
}
