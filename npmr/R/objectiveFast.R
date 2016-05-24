objectiveFast <-
function(B, P, W, lambda) {
    -sum(log(P[W])) + lambda*sum(abs(svd(B)$d))
}
