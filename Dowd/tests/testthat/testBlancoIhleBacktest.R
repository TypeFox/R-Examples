test_that("BlancoIhleBacktest.",{
  
  # Tests based on first 19 return data of SP500 returns from managers data set
  # Success - 1 
  a <- c(0.034000,  0.009300,  0.009600,  0.014700,  0.025800,  0.003800, 
         -0.044200,  0.021100,  0.056300,  0.027600,  0.075600, -0.019800,
         0.062500,  0.007800, -0.041100,  0.059700,  0.060900,  0.044800,
         0.079600) # First 19 SP500 returns from managers dataset in PerformanceAnalytics
  b <- c(0.03400, -0.01424, -0.00939, -0.00942, -0.00945, -0.00710,  0.01060,
         0.00580,  0.00100, -0.00380, -0.00435,  0.01508,  0.01272,  0.01036,
         0.03045,  0.02832,  0.02619,  0.02406,  0.02193) # Computed Using HSVaR
  # at 90% confidence level
  c <- c(0.03400, -0.00930, -0.00930, -0.00930, -0.00930, -0.00380, 0.04420,  
         0.04420,  0.04420,  0.02020,  0.02020,  0.03200, 0.03200,  0.03200,
         0.04265,  0.04265, 0.04265, 0.04265, 0.04265) # Computed using HSES at 
  # 90% confidence level
  
  expect_equal(26.5564, BlancoIhleBacktest(a, b, c, .9), tol=1)
  
})