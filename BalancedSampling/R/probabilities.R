probabilities <- function (a, n) {
  # borrowed as is from package: sampling (GPL >=2)
  nnull = length(a[a == 0])
  nneg = length(a[a < 0])
  if (nnull > 0) 
    warning("there are zero values in the initial vector a\n")
  if (nneg > 0) {
    warning("there are ", nneg, " negative value(s) shifted to zero\n")
    a[(a < 0)] = 0
  }
  if (identical(a, rep(0, length(a)))) 
    pik1 = a
  else {
    pik1 = n * a/sum(a)
    pik = pik1[pik1 > 0]
    list1 = pik1 > 0
    list = pik >= 1
    l = length(list[list == TRUE])
    if (l > 0) {
      l1 = 0
      while (l != l1) {
        x = pik[!list]
        x = x/sum(x)
        pik[!list] = (n - l) * x
        pik[list] = 1
        l1 = l
        list = (pik >= 1)
        l = length(list[list == TRUE])
      }
      pik1[list1] = pik
    }
  }
  pik1
}