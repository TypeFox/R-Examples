numericalDerivative <-
function(x, g, k, delta=.001) {
    # numerical kth derivative of g(x)
    if (k == 0) {
         g(x)
    } else {
         (numericalDerivative(x+delta, g, k-1, delta) -
             numericalDerivative(x-delta, g, k-1, delta))/
             (2*delta)
    }
}
