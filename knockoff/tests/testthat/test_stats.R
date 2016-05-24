test_that('Statistics obey antisymmetry property', {
  n = 5; p = 2;
  prob = random_problem(n, p)
  G = cbind(prob$X, knockoff.create(prob$X))
  y = prob$y
    
  i = sample(p, sample(p)) # Indices to swap.
  G_swap = G
  G_swap[,c(i,i+p)] <- G[,c(i+p,i)]
  
  expect_antisymmetric <- function(stat) {
    orig = 1:p; ko = (p+1):(2*p);
    expect_equal(stat(G[,orig],G[,ko],y), 
                 stat(G_swap[,orig],G_swap[,ko],y) * ifelse(1:p %in% i, -1, 1))
  }
  expect_antisymmetric(knockoff.stat.fs)
  expect_antisymmetric(knockoff.stat.fs_omp)
  expect_antisymmetric(knockoff.stat.lasso_difference)  
  expect_antisymmetric(knockoff.stat.lasso_signed_max)
})

test_that('Finding the max lambda in lasso works for orthonormal design', {
  n = 10; p = 5; amplitude = 3;
  X = qr.Q(qr(rnorm_matrix(n,p)))
  beta = amplitude * rnorm(p)
  y = X %*% beta + rnorm(n)
    
  beta_ls = as.vector(t(X) %*% y)
  expect_equal(lasso_max_lambda_lars(X, y), abs(beta_ls))
  expect_equal(lasso_max_lambda_glmnet(X, y, nlambda = 1e4), abs(beta_ls),
               tolerance = 1e-3)
})