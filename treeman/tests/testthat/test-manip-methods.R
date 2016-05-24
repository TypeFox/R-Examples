# LIBS
library(treeman)
library(testthat)

# TEST FUNCTIONS
randomLineage <- function(n, tree) {
  # add random monophyletic taxonyms to tree
  addLname <- function(nd, tree) {
    tree@ndlst[[nd]][['txnym']] <- lname
    pnds <- tree@ndlst[[nd]][['ptid']]
    if(!is.null(pnds)) {
      for(pnd in pnds) {
        tree <- addLname(pnd, tree)
      }
    }
    tree
  }
  nds <- tree@nds
  nds <- nds[nds != tree@root]
  lname <- paste0('l', 1)
  tree <- addLname(tree@root, tree)
  for(i in 2:(n-1)) {
    lname <- paste0('l', i)
    nd <- sample(nds, 1)
    tree <- addLname(nd, tree)
  }
  tree
}
randomTips <- function(n, tree) {
  # generate random tips with lineages for pinning
  lngs <- ends <- tip_ids <- rep(NA, n)
  nds <- names(tree@ndlst)
  nds <- nds[nds != 'n1']
  for (i in 1:n) {
    random_nd <- sample(nds, 1)
    l <- c(getNdLng(tree, random_nd),
           paste0('new_l', i))
    lngs[i] <- list(l)
    ends[i] <- runif(max=tree@age, min=0, n=1)
    tip_ids[i] <- paste0('new_', i)
  }
  list("l"=lngs, "e"=ends, "t"=tip_ids)
}

# RUNNING
context('Testing \'manip-methods\'')
# test_that('unroot() works', {
#   tree <- randTree(10)
#   tree <- unroot(tree)
#   expect_that(length(tree['root']), equals(0))
# })
test_that('addTip() works', {
  # random tree + basic stats
  data(mammals)
  tree <- getSubtree(mammals, "n821")
  pd_before <- tree['pd']
  age_before <- tree['age']
  ntips_before <- tree['ntips']
  # add random tip
  sister <- sample(tree@all[tree@all != tree@root], 1)
  sister_pd_before <- tree[[sister]]['pd']
  sister_age <- getNdAge(tree, sister)
  parent_age <- getNdAge(tree, tree@ndlst[[sister]][['prid']][1])
  start <- runif(min=sister_age, max=parent_age, n=1)
  end <- runif(min=0, max=start, n=1)
  tree <- addTip(tree, tid='new_tip', sid=sister, start=start, end=end,
                 pid='new_nd')
  tree@ndlst[['new_nd']][['spn']]
  # test if successful
  expect_that(length(tree@ndlst[[tree['root']]][['kids']]), equals(tree['ntips']))
  expect_that(validObject(tree), is_true())
  #expect_that(tree['ply'], is_false())
  expect_that(tree['age'], equals(age_before))
  expect_that(tree['ntips'], equals(ntips_before + 1))
  expect_that(tree['pd'], equals(pd_before + (start-end)))
  expect_that(tree[['new_nd']]['pd'], equals(sister_pd_before + (start - end)))
  expect_false(any(duplicated(tree[['new_nd']]['kids'])))
})
test_that('rmTip() work', {
  data(mammals)
  tree <- mammals
  pd_before <- tree['pd']
  tid <- sample(tree['tips'], 1)
  tid_spn <- getNdSlt(tree, id=tid, slt_nm='spn')
  tree <- rmTip(tree, tid)
  expect_that(tree['ntips'], equals(mammals['ntips'] - 1))
  expect_that(pd_before-tid_spn, equals(tree['pd']))
})
test_that('pinTips() work', {
  n_start <- 10
  n_add <- 20
  tree <- randTree(n_start)
  tree <- randomLineage(n_start/2, tree)
  pd_before <- tree['pd']
  age_before <- tree['age']
  rdata <- randomTips(n_add, tree)
  tree <- pinTips(tree, tids=rdata[["t"]],
                  lngs=rdata[["l"]],
                  ends=rdata[["e"]])
  expect_that(validObject(tree), is_true())
  #expect_that(tree['ntips'], equals(n_start+n_add))  # not necessarily true
  expect_that(pd_before, is_less_than(tree['pd']))
  expect_that(tree[['new_1']]['txnym'], is_a('character'))
  #expect_that(age_before, equals(tree['age']))  # not necessarily true
  writeTree(tree, file='test.tre')  # expect no error
})
if(file.exists('test.tre')) {
  file.remove('test.tre')
}