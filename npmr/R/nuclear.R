nuclear <-
function(B) {
    sum(abs(svd(B)$d))
}
