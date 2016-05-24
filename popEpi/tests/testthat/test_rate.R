context('rate')

# simultate test data
set.seed(5)
p18 <- data.table( OBS=round(runif(36)*10), PYRS=round(runif(36)*10000), AGEGROUP=1:18, COV = rep(c(1,2), each = 18))
set.seed(5)
p20 <- data.table( OBS=round(runif(20)*10), PYRS=round(runif(20)*10000), AGEGROUP=1:20, COV = rep(c(1,2), each = 20))
set.seed(5)
p101 <- data.table( OBS=round(runif(101)*10), PYRS=round(runif(101)*10000), AGEGROUP=1:101, COV = rep(c(1,2), each = 101))
p18b <- data.table(p18)
setnames(p18b, c('OBS','PYRS','AGEGROUP'), c('obs','pyrs','agegroup'))
wv <- c(.1,.1,.1,.2,.2,.2,.2,.3,.3,.4,.5,.5,.5,.4,.4,.3,.2,.1)

test_that("rate works with different weights", {
  w1 <- rate(data = p18, obs = 'OBS', pyrs = 'PYRS', print = 'COV', adjust = 'AGEGROUP', weights = c(1:18))
  w2 <- rate(data = p20, obs = 'OBS', pyrs = 'PYRS', print = 'COV', adjust = 'AGEGROUP', weights = 'world_2000_20of5')
  w3 <- rate(data = p20, obs = 'OBS', pyrs = 'PYRS', print = 'COV', adjust = 'AGEGROUP', weights = 'cohort')
  w4 <- rate(data = p20, obs = 'OBS', pyrs = 'PYRS', print = 'COV', adjust = NULL, weights = NULL)
  w5 <- rate(data = p20, obs = 'OBS', pyrs = 'PYRS', print = NULL, adjust = NULL, weights = NULL)
  
  expect_equal(sum(w1$PYRS), p18[,sum(PYRS)])
  expect_equal(sum(w2$OBS), p20[,sum(OBS)])
  expect_equal(sum(w3$OBS), p20[,sum(OBS)])
  expect_equal(w4$rate, p20[,list(sum(OBS)/sum(PYRS)), by ='COV'][, V1])
  expect_equal(w5$rate, p20[,list(sum(OBS)/sum(PYRS))][, V1])
  expect_is(w1, 'rate')
  expect_is(w2, 'rate')
  expect_is(w2, 'data.frame')
  if(getOption("popEpi.datatable")) {
    expect_is(w2, 'data.table')
  }
})

test_that("makeWeightsDT works in rate", {
  set.seed(5)
  p18 <- data.table( OBS=round(runif(36)*10), PYRS=round(runif(36)*10000), AGEGROUP=1:18, COV = rep(c(1,2), each = 18))
  op <- c('OBS',  'PYRS')
  mw1 <- makeWeightsDT(p18, adjust = substitute(factor(AGEGROUP, 1:18, 1:18)), weights = wv, print = NULL, values = list(op))
  mw2 <- makeWeightsDT(p18, adjust = NULL, weights = NULL, print = substitute(COV), values = list(op))
  
  attlist <- attr( mw1, 'makeWeightsDT')
  
  expect_equal(c(attlist$adVars, attlist$vaVars, 'weights'), c('factor','OBS','PYRS','weights') )
  
  expect_equal(mw1[,as.character(factor)], as.character(1:18))
  expect_equal(mw2[,OBS], p18[,sum(OBS), by=COV][,V1])
})


test_that("names dont cause problems", {
  w1 <- rate(data = p18b, obs = 'obs', pyrs = 'pyrs', print = 'COV', adjust = 'agegroup', weights = 'nordic')
  w2 <- rate(data = p18b, obs = 'obs', pyrs = 'pyrs', print = 'COV', adjust = 'agegroup', weights = 'cohort')
  w3 <- rate(data = p18b, obs =  obs,  pyrs =  pyrs,  print =  COV,  adjust =  agegroup,  weights = 'cohort')
  w5 <- rate(data = p18b, obs = 'obs', pyrs = 'pyrs', print = 'COV', adjust = 'agegroup', weights = 'world_1966_18of5')
  w6 <- rate(data = p20, obs = 'OBS', pyrs = 'PYRS', print = 'COV', adjust = 'AGEGROUP', weights = 'world_2000_20of5')
  w7a <- rate(data = p20, obs = OBS, pyrs = 'PYRS', print = 'COV', adjust = 'AGEGROUP', weights = 'world_2000_20of5')
  w7b <- rate(data = p20, obs = OBS, pyrs = PYRS, print = 'COV', adjust = 'AGEGROUP', weights = 'world_2000_20of5')
  
  wr <- p18b[,list(obs=sum(obs),pyrs =sum(pyrs)), by ='COV']
  
  expect_equal(w2$obs, w1$obs)
  expect_equal(w2$pyrs, w1$pyrs)
  expect_equal(wr$obs, w1$obs)
  expect_equal(wr$pyrs, w1$pyrs)
  expect_equal(w2, w3, check.attributes = FALSE)
})


test_that("rate works with different weights an subset", {
  s0 <- rate(data = p18, obs = 'OBS', pyrs = 'PYRS', print = NULL, adjust = 'AGEGROUP', weights = c(1:18))
  s0 <- rate(data = p18, obs = 'OBS', pyrs = 'PYRS', print = NULL, adjust = 'AGEGROUP', weights = c(1:18), subset = COV==1)
  s1 <- rate(data = p18, obs = 'OBS', pyrs = 'PYRS', print = 'COV', adjust = 'AGEGROUP', weights = c(1:18), subset = COV==1)
  s2 <- rate(data = p20, obs = 'OBS', pyrs = 'PYRS', print = 'COV', adjust = 'AGEGROUP',  weights = 'world_2000_20of5', subset = COV == 2)
  s3 <- rate(data = p20, obs = 'OBS', pyrs = 'PYRS', print = 'COV', adjust = 'AGEGROUP', weights = 'cohort', subset = COV==1)
  s4 <- rate(data = p20, obs = 'OBS', pyrs = 'PYRS', print = 'COV', adjust = NULL, weights = NULL, subset = AGEGROUP != 1)
  s5 <- rate(data = p20, obs = 'OBS', pyrs = 'PYRS', print = NULL, adjust = NULL, weights = NULL, subset = COV == 1)
  
  
  expect_equal(sum(s1$PYRS), p18[COV==1,sum(PYRS)])
  expect_equal(sum(s2$OBS), p20[COV==2,sum(OBS)])
  expect_equal(sum(s3$OBS), p20[COV==1,sum(OBS)])
  expect_equal(s4$rate, p20[AGEGROUP!= 1,list(sum(OBS)/sum(PYRS)), by ='COV'][, V1])
  expect_equal(s5$rate, p20[COV==1,list(sum(OBS)/sum(PYRS))][, V1])
  expect_is(s3, 'rate')
})

test_that("rate works with different weights and syntaxies", {
  
  
  wv <- c(.1,.1,.1,.2,.2,.2,.2,.3,.3,.4,.5,.5,.5,.4,.4,.3,.2,.1)
  
  s0 <- rate(data = p18, obs = 'OBS', pyrs = 'PYRS', print = NULL, adjust = NULL)
  s1 <- rate(data = p18, obs = 'OBS', pyrs = 'PYRS', adjust = factor(AGEGROUP, 1:18, 1:18), weights = c(.1,.1,.1,.2,.2,.2,.2,.3,.3,.4,.5,.5,.5,.4,.4,.3,.2,.1))
  s2 <- rate(data = p18, obs = 'OBS', pyrs = 'PYRS', adjust = factor(AGEGROUP, 1:18, 1:18), weights = wv)
  #s3 <- rate(data = p18, obs = 'OBS', pyrs = 'PYRS', adjust = factor(AGEGROUP, 1:18, 1:18), weights = NULL) ???
  expect_is(s1, 'rate')
  
  expect_equal(sum(s1$PYRS), p18[,sum(PYRS)])
  expect_equal(sum(s2$OBS), p18[,sum(OBS)])
  #expect_equal(sum(s3$OBS), p18[,sum(OBS)])
  
  
  # non working syntaxes
  expect_error( 
    rate(data = p18, obs = 'OBS', pyrs = 'PYRS', print = NULL, adjust = 'AGEGROUP', weights = list(1:18))
  ) # non named list
  expect_error(
    rate(data = p18, obs = 'OBS', pyrs = 'PYRS', adjust = factor(AGEGROUP, 1:18, 1:18), weights = list(c(.1,.1,.1,.2,.2,.2,.2,.3,.3,.4,.5,.5,.5,.4,.4,.3,.2,.1)))
  )
  expect_error(
    rate(data = p18, obs = 'OBS', pyrs = 'PYRS', adjust = list(factor(AGEGROUP, 1:18, 1:18), COV), weights = list( wv,  c(0.5,0.5)))
  ) # a list length of 2
  expect_error(
    rate(data = p18, obs = 'OBS', pyrs = 'PYRS', print = 'COV', adjust = c('AGEGROUP','COV'), weights = list(COV = c(.5,.5), AGEGROUP = 1:18))
  ) # duplicated names
  
  # working
  s10 <- rate(data = p18, obs = 'OBS', pyrs = 'PYRS', 
              adjust = list(agegr = factor(AGEGROUP, 1:18, 1:18)), 
              weights = list(agegr = c(.1,.1,.1,.2,.2,.2,.2,.3,.3,.4,.5,.5,.5,.4,.4,.3,.2,.1)))
  
  s11 <- rate(data = p18, obs = 'OBS', pyrs = 'PYRS', 
              adjust = list(agegr = factor(AGEGROUP, 1:18, 1:18), gender = factor(COV, 1:2, 1:2)), 
              weights = list(agegr = c(.1,.1,.1,.2,.2,.2,.2,.3,.3,.4,.5,.5,.5,.4,.4,.3,.2,.1), gender = c(1,1)))
  
  s12 <- rate(data = p18, obs = 'OBS', pyrs = 'PYRS', 
              adjust = list(agegr = factor(AGEGROUP, 1:18, 1:18)), 
              weights = list(agegr = wv))
  
  
  s13 <- rate(data = p18, obs = 'OBS', pyrs = 'PYRS', print = NULL, adjust = 'AGEGROUP', weights = list(AGEGROUP = 1:18))
  s14a <- rate(data = p18, obs = 'OBS', pyrs = 'PYRS', print = NULL, adjust = c('AGEGROUP','COV'), weights = list(AGEGROUP = 1:18, COV = c(.5,.5))) # SAMA1
  s14b <- rate(data = p18, obs = 'OBS', pyrs = 'PYRS', print = NULL, adjust = c('AGEGROUP','COV'), weights = list(COV = c(.5,.5), AGEGROUP = 1:18)) # SAMA1
  
  s16a <- rate(data = p20, obs = 'OBS', pyrs = 'PYRS', print = 'COV', adjust = 'AGEGROUP', weights = 'world_2000_20of5') #
  s16b <- rate(data = p20, obs = 'OBS', pyrs = 'PYRS', print = 'COV', adjust = list(AGEGROUP), weights = 'world_2000_20of5') #
  s16c <- rate(data = p20, obs = 'OBS', pyrs = 'PYRS', print = 'COV', adjust = AGEGROUP, weights = 'world_2000_20of5') #
  
  
  
  # Works
  s21 <- rate(data = p18, obs = 'OBS', pyrs = 'PYRS', print = 'COV', adjust = 'AGEGROUP', weights = (1:18), subset = COV==1)
  s22 <- rate(data = p20, obs = 'OBS', pyrs = 'PYRS', print = 'COV', adjust = 'AGEGROUP', weights = 'world_2000_20of5', subset = COV == 2)
  s23 <- rate(data = p20, obs = 'OBS', pyrs = 'PYRS', print = 'COV', adjust = NULL, weights = NULL, subset = AGEGROUP != 1)
  s24 <- rate(data = p20, obs = 'OBS', pyrs = 'PYRS', print = NULL, adjust = NULL, weights = NULL, subset = COV == 1)
  
  ## internal weights
  s23a <- rate(data = p20, obs = 'OBS', pyrs = 'PYRS', print = 'COV', adjust = 'AGEGROUP', weights = 'internal')
  s23b <- rate(data = p20, obs = 'OBS', pyrs = PYRS,   print = 'COV', adjust = 'AGEGROUP', weights = 'internal')
  s23c <- rate(data = p20, obs = OBS, pyrs = PYRS,   print = COV, adjust = AGEGROUP, weights = 'internal')
  s24a <- rate(data = p18, obs = 'OBS', pyrs = 'PYRS', adjust = list(AGEGROUP, COV), weights = "internal")
  s24b <- rate(data = p18, obs = 'OBS', pyrs = 'PYRS', adjust = c('AGEGROUP', 'COV'), weights = "internal")
  

})


test_that("rate works with missing values", {
  p18c <- copy(p18)
  p18c[c(1,6), PYRS := NA]
  expect_warning( rate(data = p18c, obs = 'OBS', pyrs = 'PYRS', adjust = 'AGEGROUP', weights = 1:18), "Data contains 2 NA values." )
})


test_that("warnings and stops works properly", {
  expect_error(
    rate(data = p18, obs = 'OBS', pyrs = 'PYRS', print = 'COV', adjust = 'AGEGROUP', weights = list(1:18, 2:19))
    )
  expect_error( stdr.weights(c('wold00_1','world66_5')) )
  expect_error( stdr.weights(c('wold00_20of5')) )
})



test_that("stdr.weights returns correct datasets", {
  al <- c('world_1966_18of5','europe','nordic',
          "world_2000_18of5","world_2000_101of1", 
          "world_2000_20of5")
  le <- c(18,18,18,18,101,20)
  expect_equal( stdr.weights(al[1])[,.N], le[1])
  expect_equal( stdr.weights(al[2])[,.N], le[2])
  expect_equal( stdr.weights(al[3])[,.N], le[3])
  expect_equal( stdr.weights(al[4])[,.N], le[4])
  expect_equal( stdr.weights(al[5])[,.N], le[5])
  expect_equal( stdr.weights(al[6])[,.N], le[6])
})
