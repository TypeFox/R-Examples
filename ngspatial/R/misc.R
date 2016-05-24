
#' Return an adjacency matrix for a square lattice.
#'
#' @details This function builds the adjacency matrix for the \code{m} by \code{n} square lattice.
#' @param m the number of rows in the lattice.
#' @param n the number of columns in the lattice. Defaults to \code{NULL}. If missing, the lattice is assumed to be \code{m} by \code{m}. 
#' @return A matrix \eqn{A} of 0s and 1s, where \eqn{A_{ij}} is equal to 1 iff vertices \eqn{i} and \eqn{j} are adjacent.
#' @export

adjacency.matrix = function(m, n = NULL)
{
    if (missing(n))
    {
        A = matrix(0, m^2, m^2)
        for (i in 1:m^2)
        {
            up = i - m
            down = i + m
            left = i - 1
            right = i + 1
            if (up > 0)
                A[i, up] = 1
            if (down <= m^2)
                A[i, down] = 1
            if (left %% m != 0)
                A[i, left] = 1
            if (i %% m != 0)
                A[i, right] = 1
        }
    }
    else
    {
        A = matrix(0, m * n, m * n)
        for (i in 1:(m * n))
        {
            up = i - n
            down = i + n
            left = i - 1
            right = i + 1
            if (up > 0)
                A[i, up] = 1
            if (down <= (m * n))
                A[i, down] = 1
            if (left %% n != 0)
                A[i, left] = 1
            if (i %% n != 0)
                A[i, right] = 1
        }
    }
    A
}

hpd = function(x, alpha = 0.05)
{
    n = length(x)
    m = round(n * alpha)
    x = sort(x)
    y = x[(n - m + 1):n] - x[1:m]
    z = min(y)
    k = which(y == z)[1]
    c(x[k], x[n - m + k])
}

is.wholenumber = function(x, tol = .Machine$double.eps^0.5)
{
    abs(x - round(x)) < tol
}

is.zero = function(x, tol = .Machine$double.eps^0.5)
{
    abs(x) < tol
}

