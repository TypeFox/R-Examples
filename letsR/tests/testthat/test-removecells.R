context("Test for .removeCells")

test_that(".removeCells works fine", {
  
  m_test <- rbind(c(1, 5, 0, 0), c(0, 5, 0, 5))
  
  expect_true(all(dim(.removeCells(m_test)) == c(1, 4)))
  
  expect_error(.removeCells(m_test[, -4]))
  
  expect_true(all(dim(.removeCells(m_test[-1, , drop=FALSE])) == c(1, 4)))
  
  
})
