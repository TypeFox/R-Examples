`strucprep` <-
function(x)
{
  distvec <- as.vector(x[lower.tri(x)])
  n <- dim(x)[1]
  dissim <- structure(distvec, Size = n, call = quote(as.dist.default(m=b)), class = "dist", Diag = FALSE, Upper = FALSE)
}

