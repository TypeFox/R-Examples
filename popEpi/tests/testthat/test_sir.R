context("SIR")


test_that("SIR w/ coh=ref=popEpi::sire", {
  ## don't skip on CRAN
  sire2 <- copy(popEpi::sire)
  sire2[, agegroup := cut(dg_age, breaks = c(0:17*5, Inf))]
  levels(sire2$agegroup) <- 1:18
  sire2[, ex_y := year(ex_date)]
  sire2[, dur := as.integer(ex_date-dg_date)/365.242199]
  ltre <- ltable(sire2[status != 0], c("ex_y", "agegroup"), 
                 expr = list(obs=.N, pyrs=sum(dur)))
  setDT(ltre)
  ltre[is.na(pyrs), pyrs := 0]
  
  sibr2 <- copy(popEpi::sibr)
  sibr2[, agegroup := cut(dg_age, breaks = c(0:17*5, Inf))]
  levels(sibr2$agegroup) <- 1:18
  sibr2[, ex_y := year(ex_date)]
  sibr2[, dur := as.integer(ex_date-dg_date)/365.242199]
  ltbr <- ltable(sibr2[status != 0], c("ex_y", "agegroup"), 
                 expr = list(obs=.N, pyrs=sum(dur)))
  setDT(ltbr)
  ltbr[is.na(pyrs), pyrs := 0]


  suppressMessages(
    sl <- sir(coh.data=ltre, coh.obs="obs", coh.pyrs="pyrs",
              ref.data=ltre, ref.obs="obs", ref.pyrs="pyrs", 
              adjust= c("agegroup","ex_y"))
  )
  ## SIR w/ coh=ref=popEpi::sire
  ## don't skip on CRAN
  expect_equal(sl$total$sir, 1)
  expect_equal(sl$total$pyrs, 13783.81, tolerance=0.01)
  expect_equal(sl$total$expected, 4595)
  expect_equal(sl$total$observed, 4595)
  
  suppressMessages(
    sl <- sir(coh.data=ltre, coh.obs="obs", coh.pyrs="pyrs",
              ref.data=ltbr, ref.obs="obs", ref.pyrs="pyrs",
              adjust= c("agegroup","ex_y"))
  )
  ## SIR w/ coh=ref=popEpi::sire"
  expect_equal(sl$total$sir, 1.39, tolerance=0.01)
  expect_equal(sl$total$pyrs, 13783.81, tolerance=0.01)
  expect_equal(sl$total$expected, 3305.04, tolerance=0.01)
  expect_equal(sl$total$observed, 4595)
})



# SIR mstate, subset + lexpand aggre ---------------------------------------

# same model
c <- lexpand( popEpi::sire[dg_date<ex_date,], status = status, birth = bi_date, exit = ex_date, entry = dg_date,
              breaks = list(per = 1990:2013, age = 0:100, fot = c(0,10,20,Inf)), 
              aggre = list(fot, agegroup = age, year = per, sex) )
# different models
c2 <- lexpand( popEpi::sire[dg_date<ex_date,], status = status, birth = bi_date, exit = ex_date, entry = dg_date,
               breaks = list(per = 1990:2010, age = 0:100, fot = c(0,10,20,Inf)), 
               aggre = list(fot, agegroup = age, year = per, sex) )



test_that("SIR works with multistate aggregated lexpand data", {
  ## don't skip on CRAN
  suppressMessages(
    suppressWarnings(  
      se <- sir( coh.data = c, coh.obs = c('from0to1','from0to2'), coh.pyrs = 'pyrs', 
                 subset = year %in% 1990:2009,
                 ref.data = data.table(popEpi::popmort), ref.rate = 'haz', 
                 adjust = c('agegroup','year','sex'), print =c('cause','fot'), mstate = 'cause')
    )
  )
  dummy.coh <- data.table(agegroup = 1:18, sex = 1, from0to1 = round(runif(18)), 
                          from0to2 = round(runif(18)+0.4),
                          pyrs = (rnorm(18)+100)*200)
  dummy.ref <- data.table(agegroup = rep(1:18,times=2), sex = 1, obs = floor(runif(36)*100),
                          categ2 = rep(c(1,2),each=18),pyrs = (rnorm(36)+100)*2000)
  suppressMessages(
    sm <- sir( coh.data = dummy.coh, coh.obs = c('from0to1','from0to2'), coh.pyrs = 'pyrs', 
               ref.data = dummy.ref, ref.obs = 'obs', ref.pyrs = 'pyrs', mstate = 'categ2',
               adjust = c('agegroup','sex','categ2'), print =c('categ2'))
  )
  # please finish this
  suppressMessages(
    s1 <- sir( coh.data = c, coh.obs = c('from0to1'), coh.pyrs = 'pyrs', subset = year %in% 1990:2009,
               ref.data = data.table(popEpi::popmort), ref.rate = 'haz', 
               adjust = c('agegroup','year','sex'), print =c('fot'))
  )
  suppressMessages(
    s2 <- sir( coh.data = c2, coh.obs = c('from0to2'), coh.pyrs = 'pyrs',
               ref.data = data.table(popEpi::popmort), ref.rate = 'haz', 
               adjust = c('agegroup','year','sex'), print =c('fot'))
  )
  s12 <- rbind( cbind(cause=1L, s1[[3]]), 
                cbind(cause=2L, s2[[3]]))
  
  # compare by hand calculated sir
  
  r <- merge(c2, data.table(popEpi::popmort), c('sex','agegroup','year'), all.x=TRUE)
  setDT(r)
  r[, exp := haz*pyrs]
  est <- r[, list(observed=sum(from0to1, na.rm=TRUE),expected=sum(exp, na.rm=TRUE)), by=.(fot)]
  est <- round(est,2)
  
  expect_is(object = se, class = 'sir')
  expect_equivalent(se[[3]], s12)
  setDT(s1[[3]])
  expect_equivalent(s1[[3]][,1:3, with=FALSE], est)
})



# SIR utils ---------------------------------------------------------------

test_that('Util functions work', {
  expect_equal( poisson.ci(5,5)$rate, 1)
})



# sir splines -------------------------------------------------------------

#library(reshape2)

test_that("SIR spline throws errors correctly", {
  skip_on_cran()
  library(splines)
  
  
  sp0 <- suppressWarnings(try(sirspline( coh.data = c, coh.obs = 'from0to2', coh.pyrs = 'pyrs',
                                         subset = year %in% 1990:2008,
                                         ref.data = data.table(popEpi::popmort), ref.rate = 'haz', 
                                         adjust = c('agegroup','year','sex'), print = NULL,
                                         spline=c('agegroup','year','fot') )))
  sp1 <- suppressWarnings(try(sirspline( coh.data = c, coh.obs = c('from0to1','from0to2'), coh.pyrs = 'pyrs',
                                         subset = year %in% 1990:2008,
                                         ref.data = data.table(popEpi::popmort), ref.rate = 'haz', 
                                         adjust = c('agegroup','year','sex'), print =c('cause'),
                                         mstate = 'cause', spline=c('agegroup','year','fot') )))
  sp2 <- suppressWarnings(try(sirspline( coh.data = c, coh.obs = c('from0to1','from0to2'), coh.pyrs = 'pyrs',
                                         subset = year %in% 1990:2008,
                                         ref.data = data.table(popEpi::popmort), ref.rate = 'haz', dependent.spline=FALSE,
                                         adjust = c('agegroup','year','sex'), print =c('cause'), 
                                         mstate = 'cause', spline=c('agegroup','year','fot') )))
  sp3 <- suppressWarnings(try(sirspline( coh.data = c, coh.obs = c('from0to1','from0to2'), coh.pyrs = 'pyrs',
                                         subset = year %in% 1990:2008,
                                         ref.data = data.table(popEpi::popmort), ref.rate = 'haz', 
                                         adjust = c('agegroup','year','sex'), print =c('cause'), 
                                         mstate = 'cause', spline='agegroup') ))
  
  sp4 <- suppressWarnings(try(sirspline( coh.data = c, coh.obs = c('from0to1','from0to2'), coh.pyrs = 'pyrs',
                                         subset = year %in% 1990:2008,
                                         ref.data = data.table(popEpi::popmort), ref.rate = 'haz', reference.points = c(2000,4),
                                         adjust = c('agegroup','year','sex'), print =c('cause'), 
                                         mstate = 'cause', spline=c('agegroup','year','fot') )))
  expect_is( object = sp1, class = 'sirspline')
  expect_is( object = sp2, class = 'sirspline')
  expect_is( object = sp3, class = 'sirspline')
  expect_is( object = sp4, class = 'sirspline')
})


test_that("print accepts a function and subset works", {
  skip_on_cran()
  suppressWarnings(
    pl1 <- sir( coh.data = c, coh.obs = c('from0to1','from0to2'), coh.pyrs = 'pyrs',
                subset = year %in% 1990:2008,
                ref.data = data.table(popEpi::popmort), ref.rate = 'haz', 
                adjust = c('agegroup','year','sex'), print = list(year.int = findInterval(year,c(1989,2000,2010))),
                mstate = 'cause')
  )
  suppressWarnings(
    pl2 <- sir( coh.data = c, coh.obs = c('from0to1','from0to2'), coh.pyrs = 'pyrs',
                subset = year %in% 1990:2008,
                ref.data = data.table(popEpi::popmort), ref.rate = 'haz', 
                adjust = c('agegroup','year','sex'), print = list(year.int = findInterval(year,c(1989,2000,2010)), sex),
                mstate = 'cause')
  )
  setDT(pl1[[2]])
  expect_equal( pl1[[2]][,year.int], 1:2)
  expect_is(object=pl2, 'sir')
})


