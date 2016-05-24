 
 
library(dplyr)


# Data generator for 2-dim anticor data
gen_data <- function(N, cor, dim = 2) {  
  rndvals <- matrix(runif(dim * N), N, dim)
  corvals <- matrix(runif(dim * N), N, dim)
  if (cor >= 0) { 
    corvals[,2:dim] <- corvals[,1]
  } else {
    corvals <- -log(corvals)
    corvals <- corvals / rowSums(corvals)
  }
  data <- as.data.frame((1 - abs(cor)) * rndvals + abs(cor) * corvals)
  names(data) <- paste0('x', 1:dim)
  return(data)
}
 


 
# Run all tests for parallel AND for non-parallel mode
for (parallelity in c(FALSE, TRUE)) {
  
  options(rPref.parallel = parallelity)

  # Run all tests for pareto and intersection
  for (`%op%` in list(`*`, `|`)) {
  
    test_that("Compare BNL and Scalagon", {
      # 2-dim anticorrelated set with some outliers
      df2 <- rbind(gen_data(1E6, -0.7, 2), data.frame(x1 = c(0,0,1,10,10), x2 = c(1,10,10,1,10)))
      # 2-dim grouped set
      df2g <- cbind(df2, data.frame(y = c(rep(1:20, 5E4), 1:5)))
      # 3-dim anticor set
      df3 <- gen_data(1E6, -0.6, 3)
      
      
      # Ignore different sorting (this is ok)
      # Deactivate Scalagon by setting alpha = 0
      
      # * Compare usual preference selection
      
      options(rPref.scalagon.alpha = 0)
      set1 <- sort(psel.indices(df2, low(x1) %op% low(x2)))
      set2 <- sort(psel.indices(df2g, low(x1) %op% low(x2)))
      set3 <- sort(psel.indices(df3, low(x1) %op% low(x2) * low(x3)))
      
      
      options(rPref.scalagon.alpha = 10)
      expect_equal(sort(psel.indices(df2, low(x1) %op% low(x2))), set1)
      expect_equal(sort(psel.indices(df2g, low(x1) %op% low(x2))), set2)
      expect_equal(sort(psel.indices(df3, low(x1) %op% low(x2) * low(x3))), set3)
      
      # * Compare top-k preference selection
      
      options(rPref.scalagon.alpha = 0)
      set1 <- sort(psel.indices(df2, low(x1) %op% low(x2), at_least = 500))
      set2 <- sort(psel.indices(df2g, low(x1) %op% low(x2), at_least = 20))
      set3 <- sort(psel.indices(df3, low(x1) %op% low(x2) * low(x3), top_level = 2))
      
      options(rPref.scalagon.alpha = 10)
      expect_equal(sort(psel.indices(df2, low(x1) %op% low(x2), at_least = 500)), set1)
      expect_equal(sort(psel.indices(df2g, low(x1) %op% low(x2), at_least = 20)), set2)
      expect_equal(sort(psel.indices(df3, low(x1) %op% low(x2) %op% low(x3), top_level = 2)), set3)
      
      
    })
  }
}