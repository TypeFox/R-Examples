objective <-
function(B, b, X, Y, lambda) {
    -logL(B, b, X, Y) + lambda*nuclear(B)
}
