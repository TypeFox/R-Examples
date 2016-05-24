# LIBS
library(treeman)
library(testthat)

# RUNNING
context('Testing \'calc-methods\'')
test_that('calc_Blnc() works', {
  tree <- readTree(text="((t1:1.0,t2:1.0):1.0,(t3:1.0,t4:1.0):1.0);")
  bstats <- calcNdsBlnc(tree, tree['nds'])
  expect_that(sum(bstats), equals(0))
})
test_that('calcDstTrp() works', {
  tree_1 <- randTree(10)
  # small chance that they will both be the same
  # tree_2 <- randTree(10)
  # res <- calcDstTrp(tree_1, tree_2, nrmlsd=TRUE)
  # expect_that(res, is_more_than(0))
  res <- calcDstTrp(tree_1, tree_1, nrmlsd=TRUE)
  expect_that(res, equals(0))
})
test_that('calcOvrlp() works', {
  tree <- readTree(text="((t1:1.0,t2:1.0):1.0,t3:1.0);")
  ovrlp <- calcOvrlp(tree, ids_1=tree['tips'], ids_2=c('t3'), nrmlsd=TRUE)
  expect_that(ovrlp, equals(1/4))
})
test_that('calcDstBLD() works', {
  tree_1 <- readTree(text="((t1:1.0,t2:1.0):1.0,t3:1.0);")
  tree_2 <- readTree(text="((t3:1.0,t2:1.0):1.0,t1:1.0);")
  d <- calcDstBLD(tree_1, tree_2, TRUE)
  expect_that(d, equals(1))
  d <- calcDstBLD(tree_1, tree_1, TRUE)
  expect_that(d, equals(0))
})
test_that('calcDstRF() works', {
  tree_1 <- readTree(text="((t1,t2),t3);")
  tree_2 <- readTree(text="((t3,t2),t1);")
  d <- calcDstRF(tree_1, tree_2, TRUE)
  expect_that(d, equals(1))
  d <- calcDstRF(tree_1, tree_1, TRUE)
  expect_that(d, equals(0))
})
test_that('calcPhyDv() works', {
  tree <- randTree(10)
  tips <- sample(tree['tips'], 3)
  pd <- calcPhyDv(tree, tips)
  parent <- getPrnt(tree, ids=tips)
  test_that(pd, is_less_than(tree@ndlst[[parent]][['pd']]))
  # add a tip with a specified length.
  sister <- sample(tips, 1)
  sister_age <- getNdAge(tree, sister)
  parent_age <- getNdAge(tree, tree@ndlst[[sister]][['prid']][1])
  start <- runif(min=sister_age, max=parent_age, n=1)
  end <- runif(min=0, max=start, n=1)
  tree <- addTip(tree, tid='new_tip', sid=sister, start=start, end=end)
  new_pd <- calcPhyDv(tree, c(tips, 'new_tip'))
  test_that(new_pd, equals(pd + (start - end)))
})
test_that('calcFrPrp() works', {
  tree <- randTree(10)
  ed_values <- calcFrPrp(tree, tree['tips'])
  expect_that(sum(ed_values), equals(tree['pd']))
})
test_that('calcDstMtrx() works', {
  tree <- randTree(10)
  ids <- tree['all']
  dmtrx <- calcDstMtrx(tree, ids)
  rndnd <- sample(ids, 1)
  expect_that(dmtrx[rndnd, rndnd], equals(0))
  expect_that(sum(dmtrx['n1', ] == tree['age']), is_more_than(0))
})