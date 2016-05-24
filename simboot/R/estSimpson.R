estSimpson <-
function (x) 
{
  n <- sum(x)
  estp <- x/n
  Simp <- Simpson(estp) * n/(n - 1)
  return(Simp)
}

