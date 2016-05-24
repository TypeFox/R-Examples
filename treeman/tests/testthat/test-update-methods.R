# LIBS
library(treeman)
library(testthat)

# RUNNING
context('Testing \'TreeMan Class\'')
test_that('.updateTreeSlts() works', {
  tree <- randTree(100)
  tree@ndlst[['t1']][['spn']] <- NULL
  tree <- treeman:::.updateTreeSlts(tree)
  expect_false(tree@wspn)
})
test_that('.updateKids() works', {
  tree <- randTree(10)
  kids <- getNdsKids(tree, tree['nds'])
  ndlst <- tree@ndlst
  tid <- sample(tree['tips'], 1)
  for(i in 1:length(ndlst)) {
    bool <- ndlst[[i]][['kids']] != tid
    ndlst[[i]][['kids']] <- ndlst[[i]][['kids']][bool]
  }
  ndlst <- treeman:::.updateKids(ndlst, tid=tid,
                       rid=tree['root'])
  tree@ndlst <- ndlst
  new_kids <- getNdsKids(tree, tree['nds'])
  res <- rep(NA, length(kids))
  for(i in 1:length(res)) {
    res[i] <- all(kids[[i]] %in% new_kids[[i]]) &
      length(kids[[i]]) == length(new_kids[[i]])
  }
  expect_true(all(res))
})
test_that('.dwndateKids() works', {
  tree <- randTree(10)
  ndlst <- tree@ndlst
  tid <- sample(tree['tips'], 1)
  ndlst <- treeman:::.dwndateKids(ndlst, tid, tree['root'])
  test_bool <- rep(NA, length(ndlst))
  for(i in 1:length(ndlst)) {
    test_bool[i] <- tid %in% ndlst[[i]][['kids']]
  }
  expect_false(any(test_bool))
})
test_that('.globalUpdateKids() works', {
  tree <- randTree(10)
  kids <- getNdsKids(tree, tree['nds'])
  ndlst <- tree@ndlst
  for(i in 1:length(ndlst)) {
    ndlst[[i]][['kids']] <- NULL
  }
  ndlst <- treeman:::.globalUpdateKids(ndlst)
  tree@ndlst <- ndlst
  new_kids <- getNdsKids(tree, tree['nds'])
  res <- rep(NA, length(kids))
  for(i in 1:length(res)) {
    res[i] <- all(kids[[i]] %in% new_kids[[i]]) &
      length(kids[[i]]) == length(new_kids[[i]])
  }
  expect_true(all(res))
})
test_that('.updateTip() works', {
  tree <- randTree(10)
  ndlst <- tree@ndlst
  tid <- sample(tree['tips'], 1)
  ndlst <- treeman:::.dwndateTip(ndlst, tid=tid, rid=tree['root'])
  new_ndlst <- treeman:::.updateTip(ndlst, tid=tid, rid=tree['root'])
  tree <- new('TreeMan', ndlst=new_ndlst, root=tree['root'])
  tree <- treeman:::.updateTreeSlts(tree)
  expect_that(tree['ntips'], equals(10))
})
test_that('.dwndateTip() works', {
  tree <- randTree(10)
  ndlst <- tree@ndlst
  tid <- sample(tree['tips'], 1)
  new_ndlst <- treeman:::.dwndateTip(ndlst, tid=tid, rid=tree['root'])
  pd_res <- kid_res <- rep(NA, length(ndlst))
  for(i in 1:length(ndlst)) {
    kid_res[i] <- !(tid %in% new_ndlst[[i]][['kids']])
    pd_res[i] <- ndlst[[i]][['pd']] > new_ndlst[[i]][['pd']]
  }
  expect_true(all(kid_res))
  expect_that(sum(pd_res), is_more_than(0))
})
test_that('.dwndateNd() works', {
  tree <- randTree(10)
  ndlst <- tree@ndlst
  nid <- sample(tree['nds'][tree['nds'] != tree['root']], 1)
  new_ndlst <- treeman:::.dwndateNd(ndlst, nid=nid, rid=tree['root'])
  pd_res <- kids_res <- rep(NA, length(ndlst))
  for(i in 1:length(ndlst)) {
    pd_res[i] <- ndlst[[i]][['pd']] > new_ndlst[[i]][['pd']]
    kids_res[i] <- nid %in% ndlst[[i]][['kids']]
  }
  expect_that(sum(pd_res), is_more_than(0))
  expect_that(all(kids_res), is_false())
})
test_that('.updateNd() works', {
  tree <- randTree(10)
  age_before <- tree['age']
  pd_before <- tree['pd']
  ntips_before <- tree['ntips']
  ndlst <- tree@ndlst
  nid <- sample(tree['nds'][tree['nds'] != tree['root']], 1)
  ndlst <- treeman:::.dwndateNd(ndlst, nid=nid, rid=tree['root'])
  new_ndlst <- treeman:::.updateNd(ndlst, nid=nid, rid=tree['root'])
  tree <- new('TreeMan', ndlst=new_ndlst, root=tree['root'])
  tree <- treeman:::.updateTreeSlts(tree)
  expect_that(tree['ntips'], equals(ntips_before))
  expect_that(tree['age'], equals(age_before))
  expect_that(tree['pd'], equals(pd_before))
})
test_that('.globalUpdateAll() works', {
  tree <- randTree(10)
  ndlst <- tree@ndlst
  for(i in 1:length(ndlst)) {
    ndlst[[i]][['prdst']] <- ndlst[[i]][['kids']] <- NULL
    ndlst[[i]][['pd']] <- 0
  }
  new_ndlst <- treeman:::.globalUpdateAll(ndlst=ndlst)
  tree <- new('TreeMan', ndlst=new_ndlst, root=tree['root'])
  tree <- treeman:::.updateTreeSlts(tree)
  expect_that(tree['ntips'], equals(10))
})