##' Compute external- or `cross'- product of 3D vectors.
##' 
##' Computes the external product \deqn{ }{ (x2 * y3 - x3 * y2, x3 * y1 - x1 *
##' y3, x1 * y2 - x2 * y1) }\deqn{ \left(x_2 y_3 - x_3 y_2,\; x_3 y_1 - x_1
##' y_3,\; x_1 y_2 - x_2 y_1 \right) }{ (x2 * y3 - x3 * y2, x3 * y1 - x1 * y3,
##' x1 * y2 - x2 * y1) }\deqn{ }{ (x2 * y3 - x3 * y2, x3 * y1 - x1 * y3, x1 *
##' y2 - x2 * y1) } of the 3D vectors in \bold{x} and \bold{y}.
##' 
##' 
##' @param x \code{n}-by-3 matrix. Each row is one \bold{x}-vector
##' @param y \code{n}-by-3 matrix. Each row is one \bold{y}-vector
##' @return \code{n}-by-3 matrix
##' @author Raoul Grasman
##' @keywords arith math array
##' @export
"extprod3d" <-
function (x, y) 
{
    x = matrix(x, ncol = 3)
    y = matrix(y, ncol = 3)
    drop(cbind(x[, 2] * y[, 3] - x[, 3] * y[, 2], x[, 3] * y[, 
        1] - x[, 1] * y[, 3], x[, 1] * y[, 2] - x[, 2] * y[, 
        1]))
}
