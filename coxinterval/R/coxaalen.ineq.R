### construct inequality constraint matrix (requires Matrix)
coxaalen.ineq <- function(w, k)
{
  if (!length(dim(w))) w <- matrix(w)
  m <- ncol(w)
  if (m == 0) mat <- matrix(1)
  else {
    val <- apply(w, 2, function(y) c(floor(min(y)), ceiling(max(y))))
    ind <- (sapply(0:(2^m - 1), function(y) as.logical(intToBits(y))[1:m])
            + (2^(0:(m - 1)) + c(0, rep(1, m - 1))))
    mat <- cbind(1, t(matrix(as.vector(val)[ind], nrow = m)))
  }
  dia <- as.matrix(do.call(Matrix::bdiag,
                           lapply(vector("list", k), function(y) mat)))
  dia - rbind(matrix(0, nrow = nrow(mat), ncol = ncol(dia)),
              dia[1:(nrow(dia) - nrow(mat)), ])
}
