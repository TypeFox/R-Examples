ed <- function(s){
  # Gives SORTED and NORMALIZED eigenvectors in u
  
  # The rank is determined numerically, dropping eigenvalues
  # that are less then a predetermined tolerance (1e-9).
  # Associated eigenvectors are eliminated.
  out = list()
  tol = 2.2204e-016
  tmp = eigen(s)
  d = tmp$values
  u = tmp$vectors
  i = abs(d)>tol
  u = u[,i]
  d = d[i]
  p = length(d)
  tmp = sort(d,index.return=TRUE)
  d =tmp$x
  i = tmp$ix
  j=p:1
  u = u[,i[j]]
  u = nrm(u) 
  out$u = u
  out
}

