library(coop)
m <- 30
n <- 10


test <- function(dense, sparse)
{
  t1 <- coop::cosine(dense)
  t2 <- coop::cosine(sparse)
  
  stopifnot(all.equal(t1, t2))
}


if (require(slam))
{
  library(slam)
  set.seed(1234)
  
  ### Very sparse, has column of 0's
  x <- coop:::dense_stored_sparse_mat(m, n, prop=.05)
  coo <- as.simple_triplet_matrix(x)
  test(x, coo)
  
  colnames(x) <- sample(letters, size=n, replace=TRUE)
  coo <- as.simple_triplet_matrix(x)
  test(x, coo)
  
  ### Not very sparse
  x <- coop:::dense_stored_sparse_mat(m, n, prop=.25)
  coo <- as.simple_triplet_matrix(x)
  test(x, coo)
  
  colnames(x) <- sample(letters, size=n, replace=TRUE)
  coo <- as.simple_triplet_matrix(x)
  test(x, coo)
}



if (require(Matrix))
{
  library(Matrix)
  set.seed(1234)
  
  ### Very sparse, has column of 0's
  x <- coop:::dense_stored_sparse_mat(m, n, prop=.05)
  csc <- as(x, "sparseMatrix")
  test(x, csc)
  
  colnames(x) <- sample(letters, size=n, replace=TRUE)
  csc <- as(x, "sparseMatrix")
  test(x, csc)
  
  ### Not very sparse
  x <- coop:::dense_stored_sparse_mat(m, n, prop=.25)
  csc <- as(x, "sparseMatrix")
  test(x, csc)
  
  colnames(x) <- sample(letters, size=n, replace=TRUE)
  csc <- as(x, "sparseMatrix")
  test(x, csc)
}
