pairwise.product <-
function(x)
{
  p <- ncol(x)
  product <- combn(seq_len(p), 2L, function(i) x[, i[1L]] * x[, i[2L]])
  product
}
