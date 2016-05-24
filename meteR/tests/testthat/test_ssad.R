context('test of ssad predictions')
## values from newman et al. 

vionut <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2,
            2, 2, 2, 3, 3, 3, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 7, 9)

ourSSF <- meteSSF(n0=101, A=16/64, A0=16)
ourSSAD <- ssad(ourSSF)

test_that('class from SSF is correct', {
  expect_is(ourSSF, 'meteSSF')
})

test_that('class from ssad is correct', {
  expect_is(ourSSAD, 'ssad')
})

test_that('predicted rank ssad is correct', {
  expect_equal(sort(meteDist2Rank(ourSSAD)), vionut)
})
