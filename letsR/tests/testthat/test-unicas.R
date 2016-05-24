context("Test for .unicas")


m_test <- matrix(runif(100), ncol=10)
colnames(m_test) <- c("a", "a", "b", "b", "b", 
                      "c", "d", "Bruno", "a1", 
                      "A")


test_that(".unicas works fine", {
  
  resu_test <- .unicas(m_test)
  expect_true(ncol(resu_test) == 7)
  expect_true(is.matrix(resu_test))
  
  
})
