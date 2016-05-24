context('exports')

tst.sum <- lrequire('single-function-encapsulation')

test_that('exposed function is callable and functions correctly', {
  expect_equal(tst.sum(1, 2), 3)
})
