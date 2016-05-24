Perms <-
function(n)#generated matrix of all possible combos of predictors in tree
{
  mat<-matrix(0, nrow=2^n, ncol=n)
  for(i in 1:n)
    {
    mat[,i]<-rep(0:1, times=2^(i-1), each=2^(n-i))
    }
  mat
}
