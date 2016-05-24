AlphaCross <-
function (L1, R1, n1, est2, k2)
{
  alpha = matrix(0, n1, k2, dimnames=list(1:n1, row.names(est2$alpha)))
  for (i in 1:n1)
  {
    tempint = L1[i] < est2$intmap[2,] & R1[i] >= est2$intmap[2,]
    alpha[i, tempint] = 1
  }
  alpha
}