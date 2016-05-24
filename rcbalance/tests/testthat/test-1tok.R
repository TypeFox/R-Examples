library(rcbalance)
library(optmatch)
context('1-to-k matching')

data(nuclearplants)
my.dist.struct <- build.dist.struct(z = nuclearplants$pr, 
	X = subset(nuclearplants[c('date','t1','t2','cap','bw','cum.n')]),
	calip.option = 'none')



test_that("1-to-k match produces k columns in matches matrix", {
  my.out <- rcbalance(my.dist.struct, k =2)
  expect_equal(ncol(my.out$matches), 2)
})

test_that("1-to-k match with refined balance has correct column count and table sums ", {
  nuke.ctrl <- nuclearplants[nuclearplants$pr ==0,]
  nuke.treat <- nuclearplants[nuclearplants$pr ==1,]
  my.out <-rcbalance(my.dist.struct, k =2, treated.info = nuke.treat, control.info = nuke.ctrl, fb.list= list('ct',c('ct','bw'))) 
  expect_equal(ncol(my.out$matches), 2)
  #count people in fb.tables, make sure ratio is right
  tableSums <- t(sapply(my.out$fb.tables, colSums))
  expect_equal(tableSums[,1],2*tableSums[,2])
})

