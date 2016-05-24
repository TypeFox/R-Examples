# LIBS
library(treeman)
library(testthat)

# RUNNING
context('Testing \'read-write-methods\'')
test_that('readTree([w/ spans]) works', {
  text <- "(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5);"
  tree <- readTree(text=text)
  expect_that(tree['pd'], equals(0.1+0.2+0.3+0.4+0.5))
  expect_that(tree['ntips'], equals(4))
  expect_that(tree[['C']]['prdst'], equals(0.3+0.5))
})
test_that('readTree([w/o spans]) works', {
  text <- "(A,B,(C,D)E);"
  tree <- readTree(text=text)
  expect_that(tree['ntips'], equals(4))
  expect_false(tree['wspn'])
})
test_that('writeTree() works', {
  t1 <- randTree(100)
  writeTree(t1, 'test.tre')
  t2 <- readTree('test.tre')
  expect_that(t1['ntips'], equals(t2['ntips']))
  expect_that(t1['nnds'], equals(t2['nnds']))
  expect_that(t1['pd'], equals(t2['pd']))
  expect_that(t1['age'], equals(t2['age']))
})
if(file.exists('test.tre')) {
  file.remove('test.tre')
}