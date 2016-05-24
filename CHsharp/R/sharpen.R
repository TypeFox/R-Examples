sharpen <-
function(x, y, lambda, B) {
    # result of sharpening operation using penalty summarized by matrix B
    # x, y are vectors of n original observations
    # output is the vector of sharpened y values
    solve(diag(rep(1, length(y)))+lambda*B%*%t(B), y)
}
