# LIBS
library(treeman)
library(testthat)

# RUNNING
context('Testing \'get-methods\'')
test_that('get_Sstr() works', {
  tree <- randTree(10)
  tips <- tree['tips']
  fwrd <- getNdsSstr(tree, tips)
  rvrse <- getNdsSstr(tree, fwrd)
  expect_that(rvrse, equals(tips))
})
test_that('getTxnyms() works', {
  data('mammals')
  nid <- sample(mammals['nds'], 1)
  txnym <- mammals[[nid]]['txnym']
  res <- getTxnyms(mammals, txnym)
  expect_true(nid %in% res[[1]])
})
test_that('getOtgrp() works', {
  tree <- randTree(10)
  rnd_nd <- sample(tree['nds'][tree['nds'] != tree['root']], 1)
  ingrp <- getNdKids(tree, rnd_nd)
  otgrp <- sample(tree['tips'][!tree['tips'] %in% ingrp], 1)
  res <- getOtgrp(tree, ids=c(ingrp, otgrp))
  expect_that(res, equals(otgrp))
})
test_that('get_Slt() works', {
  tree <- randTree(10)
  nd_spns <- getNdsSlt(tree, slt_nm="spn",
                       ids=tree['all'])
  expect_that(sum(nd_spns), equals(tree['pd']))
})
test_that('get_Kids() works', {
  tree <- randTree(10)
  kids <- getNdsKids(tree, tree['nds'])
  expect_true(all(kids$n1 %in% paste0("t", 1:10)))
})
test_that('get_Age() works', {
  tree <- randTree(10)
  root_age <- tree['age']
  nd_ages <- getNdsAge(tree, tree['nds'])
  expect_true(all(nd_ages <= root_age))
})
test_that('getSpnAge() works', {
  tree <- randTree(10)
  tip_age <- getSpnAge(tree, sample(tree['tips'], 1))
  expect_that(tip_age[['start']], is_more_than(tip_age[['end']]))
})
test_that('getSpnsAge() works', {
  tree <- randTree(10)
  tip_ages <- getSpnsAge(tree, tree['tips'])
  res <- all(tip_ages[ ,'start'] > tip_ages[ ,'end'])
  expect_true(res)
})
test_that("getPrnt() works", {
  tree <- readTree(text="(((A,B),(C,D)),(E,F));")
  prnt <- getPrnt(tree, ids=c("A", "C"))
  expect_true(prnt == "n2")
})
test_that("getPath() works", {
  tree <- randTree(10)
  pth <- getPath(tree, from="t1", to="t10")
  prnt <- getPrnt(tree, ids=c('t1', "t10"))
  expect_true(prnt %in% pth)
  expect_that(pth[1], equals('t1'))
  expect_that(pth[length(pth)], equals('t10'))
})
test_that("get_Prid() works", {
  tree <- randTree(10)
  prid <- getNdPrid(tree, id='n1')
  expect_that(prid, is_null())
  prids <- getNdsPrid(tree, tree['nds'])
  lst_nds <- unlist(lapply(prids, function(n) n[length(n)]))
  expect_true(all(lst_nds == "n1"))
})
test_that("get_Ptid() works", {
  tree <- randTree(10)
  ptids <- getNdsPtid(tree, tree['nds'])
  n1_ptids <- tree['all'][tree['all'] != 'n1']
  expect_true(all(n1_ptids %in% ptids[['n1']]))
  expect_that(ptids[['t1']], is_null())
})
test_that("get_Lng() works", {
  # TODO: redo with setNds
  tree <- randTree(10)
  for(i in 1:length(tree@ndlst)) {
    tree@ndlst[[i]]$txnym <- paste0("l", sample(1:1000, 1))
  }
  lngs <- getNdsLng(tree, tree['tips'])
  rnd1 <- sample(1:length(lngs), 1)
  rnd2 <- sample(1:length(lngs), 1)
  expect_that(sum(lngs[[rnd1]] %in% lngs[[rnd2]]),
              is_more_than(0))
})
test_that("getSubtree() works", {
  tree <- randTree(10)
  subtree <- getSubtree(tree, 'n2')
  expect_that(tree['ntips'], is_more_than(subtree['ntips']))
  expect_that(tree['nnds'], is_more_than(subtree['nnds']))
  expect_that(tree['pd'], is_more_than(subtree['pd']))
  expect_that(tree['age'], is_more_than(subtree['age']))
})