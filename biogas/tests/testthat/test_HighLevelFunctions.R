# test values form the high level functions
# Charlotte
# 22 July 2015
# modified 14 April 2016 SDH

context("Tests high level functions")

# predBg
test_that('default methane prediction from COD is stable', {
  expect_equal(round(predBg(COD = 1),4), 349.3842)
})

test_that('default methane prediction from COD is stable with fs argument', {
  expect_equal(signif(predBg(COD = 1, fs = 0.1), 6), signif(349.3842*0.9, 6))
})

test_that('default methane prediction from COD is stable with fd argument', {
  expect_equal(signif(predBg(COD = 1, fd = 0.9), 6), signif(349.3842*0.9, 6))
})

test_that('default methane prediction from COD is stable', {
  expect_equal(round(predBg(COD = 1),4), 349.3842)
})

test_that('default methane prediction from predBg using a formula is stable', {
  expect_equal(predBg('C6H10O6'), calcCOD('C6H10O6')*predBg(COD = 1) )
})


# cumBg.R
# NTS: for the moment everything is tested once, it could be good to separatae tests for - cH4 volume vs biogas volume and rates
# NTS: might also be good to store the test data framse somewhere

# volume molar values:
vmch4 <- 22360.588
vmco2 <- 22263.009

# cumulative sum and rates
# volumetric
test_that("cumulative sum and rates are corrrectly calculated with default values", {
  test.vol <- data.frame(id = rep(paste0('R_',1),5), time= c(2, 4, 5, 1, 3), vol = rep(20, 5))
  res <- data.frame(id = rep('R_1', 6), time = c(0:5), vol = c(0, rep(20, 5)), vBg = c(0, rep(20,5)), cvBg = c(0, 20, 40, 60, 80, 100), rvBg = c(NA,rep(20,5)))

  expect_equal( cumBg(dat = test.vol), res)
})

test_that("cumulative sum and rates are corrrectly calculated with default values , except addt0 = FALSE", {
  test.vol <- data.frame(id = rep(paste0('R_',1),5), time= c(2, 4, 5, 1, 3), vol = rep(20, 5))
  res <- data.frame(id = rep('R_1', 5), time = c(1:5), vol = rep(20, 5), vBg = rep(20,5), cvBg = c(20, 40, 60, 80, 100), rvBg = c(NA,rep(20,4)))
  
  expect_equal( cumBg(dat = test.vol, addt0 = FALSE), res)
})


test_that("cumulative sum, rates and methane volume are corrrectly calculated using one value for comp", {
  test.vol <- data.frame(id = rep(paste0('R_',1),5), time= c(2, 4, 5, 1, 3), vol = rep(20, 5))
# wanted result : calculation for each column from data in test.mass
  res <- data.frame(id = rep('R_1', 6), time = c(0:5), vol = c(0, rep(20, 5)), xCH4 = c(NA, rep(0.6, 5)))
  res$vBg <- stdVol(c(0, rep(20,5)), temp = 35, pres = 1)
  res$vCH4 <- res$vBg*0.6*vmch4 / (0.6*vmch4 + 0.4*vmco2) 
  res$cvBg <- rep(0,6)
 for(i in 2:6){
  res[i, 'cvBg'] = res[i, 'vBg'] + res[i-1, 'cvBg'] 
}
  res$cvCH4 <- res$cvBg*0.6*vmch4/ (0.6*vmch4 + 0.4*vmco2)
  res$rvBg <- c(NA, res[2:6, 'vBg']) # because time interval = 1
  res$rvCH4 <- c(NA, res[2:6, 'vCH4'])
  
  expect_equal( a <- cumBg(dat = test.vol, comp = 0.6, temp = 35, pres = 1 ), res)
})


# gravimetric
test_that("cumulative sum, rates and methane volume are corrrectly calculated using mass + one required value for comp, temp and pres", {
  # set up data frame to test the function
  test.mass <- data.frame(id = rep(paste0('R_',1),6), time= c(2, 4, 5, 1, 3, 0), mass = c(95, 85, 80, 100, 90, 105))
  # wanted result : calculation for each column from data in test.mass
  res <- data.frame(id = rep('R_1',6) , time = c(0:5), mass = c(105, 100, 95, 90, 85, 80), xCH4 = rep(0.6, 6),  
                    massloss = c(0,rep(5,5)), cmassloss = c(0, 5, 10, 15, 20, 25))
  res <- data.frame( res, mass2vol(res$massloss, xCH4 = 0.6, temp = 35 , pres = 1, value= 'all'))
  res <- res[,-match('vCO2',names(res))]
  res$cvBg <- rep(0,6)
  for(i in 2:6){
  res[i, 'cvBg'] = res[i, 'vBg'] + res[i-1, 'cvBg'] 
  }
  res$cvCH4 <- res$cvBg*0.6*vmch4/ (0.6*vmch4 + 0.4*vmco2)
  res$rvBg <- c(NA, res[2:6, 'vBg']) # because time interval = 1
  res$rvCH4 <- c(NA, res[2:6, 'vCH4'])
  
  expect_equal( cumBg(dat = test.mass, dat.type = "mass", comp = 0.6, temp = 35 , pres = 1), res)
})

# NTS : maybe change order of thest in grav : or test all at the same time?
# NTS check message content in more details
test_that("missing arguments in gravimetric throw an error", {
  test.mass <- data.frame(id = rep(paste0('R_',1),6), time= c(2, 4, 5, 1, 3, 0), mass = c(95, 85, 80, 100, 90, 105))
  expect_that(cumBg(dat = test.mass, dat.type = "mass"), throws_error())
  expect_that(cumBg(dat = test.mass, dat.type = "mass", temp = 35), throws_error())
  expect_that(cumBg(dat = test.mass, dat.type = "mass", temp = 35, pres =1), throws_error())
})

test_that("cumBg gives a message error if cumulative mass is negative", {
  test.mass <- read.csv('test.mass.csv')
  expect_that(cumBg(dat = test.mass, dat.type = "mass", comp = "comp", temp = 35, pres = 1, time.name = 'days'), throws_error())
  
}
)


# summBg
# uses same res as volumetric testing n°2
test_that("summBg correctly calculates means", {
  # same result data frame as cumBg
  res <- data.frame(id = rep('R_1', 6), time = c(0:5), vol = c(0, rep(20, 5)), xCH4 = c(NA, rep(0.6, 5)))
  res$vBg <- stdVol(c(0, rep(20,5)), temp = 35, pres = 1)
  res$vCH4 <- res$vBg*0.6*vmch4 / (0.6*vmch4 + 0.4*vmco2) 
  res$cvBg <- rep(0,6)
  for(i in 2:6){
    res[i, 'cvBg'] = res[i, 'vBg'] + res[i-1, 'cvBg'] 
  }
  res$cvCH4 <- res$cvBg*0.6*vmch4/ (0.6*vmch4 + 0.4*vmco2)
  res$rvBg <- c(NA, res[2:6, 'vBg']) # because time interval = 1
  res$rvCH4 <- c(NA, res[2:6, 'vCH4'])
  
  # add a reactor to res data frame 
  test.sum1 <- cbind(res[, 1:4], res[, -c(1:4)]*1.2)
  test.sum1$vol <- test.sum1$vol*1.2
  test.sum1$id <- rep('R_2', 6)
  
  test.sum2 <- cbind(res[, 1:4], res[, -c(1:4)]*3)
  test.sum2$vol <- test.sum2$vol*3
  test.sum2$id <- rep('R_3', 6)
  
  test.sum3 <- cbind(res[, 1:4], res[, -c(1:4)]*3.3)
  test.sum3$vol <- test.sum3$vol*3.3
  test.sum3$id <- rep('R_4', 6)
  
  test.sum <- rbind(res, test.sum1)
  test.sum <- rbind(test.sum, test.sum2)
  test.sum <- rbind(test.sum, test.sum3)

  setup <- data.frame( id = c('R_1', 'R_2', 'R_3', 'R_4'), descrip = c('inoc', 'inoc', 'A', 'A')) 
  
  res.test.sum <- 
    data.frame (descrip = c('inoc', 'A'), 
                time = c(5,5), 
                mean = 
                  c(mean(test.sum[test.sum$time == 5 & test.sum$id %in% c('R_1', 'R_2'), 'cvCH4']), 
                    mean(test.sum[test.sum$time == 5 & test.sum$id %in% c('R_3', 'R_4'), 'cvCH4']))
                )
  expect_equal(summBg(test.sum, setup, when = 5, sort = FALSE)[, 1:3], res.test.sum)

}) 


