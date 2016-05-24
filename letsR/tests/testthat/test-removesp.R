context("Test for .removeSp")

test_that(".removeSp works fine", {
  
  m_test <- rbind(c(1, 5, 0, 5), c(0, 5, 0, 5))
  
  expect_true(all(dim(.removeSp(m_test)) == c(2, 3)))
  
  expect_error(.removeSp(m_test[, -4]))
  
  expect_true(all(dim(.removeSp(m_test[-1, , drop=FALSE])) == c(1, 3)))
  
})

