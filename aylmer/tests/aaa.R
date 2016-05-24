require(aylmer)
options(warn = -1)

f <-
  function(x,tol=1e-6)
{
  jj1 <- fisher.test(x)$p.value
  jj2 <- aylmer.test(x)$p.value
  stopifnot(abs( (jj1-jj2)/(jj1+jj2)) < tol)
}

data(iqd)

# Following two tests inspired by an anonymous JSS referee:
f(shifts[-1,]) 
f(shifts[,-3])

# Following tests made up by me:
f(diag(3 + 1:3))
f(matrix(1,3,4))
