pairwise_association <- function(mat, method="condentropy")
{
  nsp <- nrow(mat)
  m_spxsp <- outer(1:nsp, 1:nsp, FUN=Vectorize(function(i,j) { do.call(method, args=list(X=mat[i,], Y=mat[j,])) }))
  
  diag(m_spxsp) <- NA
  
  return(m_spxsp)
}