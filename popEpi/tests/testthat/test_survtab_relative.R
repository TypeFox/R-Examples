context("popEpi::survtab vs. relsurv::rs.surv")

# survtab vs. relsurv::rs.surv --------------------------------------------
test_that("relative survivals about the same as relsurv's", {
  
  library(survival)
  library(relsurv)
  library(Epi)
  
  # male
  pm <- data.table(popEpi::popmort)
  # pm[, surv := 1L]
  pm[, surv := exp(-haz)]
  pm.m <- cast_simple(pm[sex==0], columns = 'year', rows = 'agegroup',  values='surv')
  pm.m[,agegroup := NULL]
  pm.m <- as.matrix(pm.m)
  # female
  pm.f <- cast_simple(pm[sex==1], columns = 'year', rows = 'agegroup',  values='surv')
  pm.f[,agegroup := NULL]
  pm.f <- as.matrix(pm.f)
  
  popm <- transrate(pm.m, pm.f, yearlim = c(1951, 2013), int.length = 1)
  
  pm[, surv := NULL]
  
  sire2 <- sire[dg_date<ex_date, ]
  sire2[, Tstop  := as.integer(ex_date - dg_date)]
  sire2[, dg_age := as.integer(dg_date - bi_date)]
  
  x <- Lexis(entry = list(age = dg_age, per = dg_date, fot = 0L),
             exit = list(fot = Tstop), 
             exit.status = as.integer(status %in% 1:2),
             entry.status = 0L, data = sire2)
  setDT(x)
  setattr(x, "class", c("Lexis","data.table", "data.frame"))
  
  ## rs.surv
  ## sex must be coded c(1,2) (male, female)
  x[, sex := 2L]
  rs.e2 <- rs.surv(Surv(lex.dur, lex.Xst!=0) ~ 1 + ratetable(age=age, sex=sex, year=per),
                   ratetable = popm, data = x, method = 'ederer2', type = "fleming-harrington", fin.date=ex_date)
  rs.pp <- rs.surv(Surv(lex.dur, lex.Xst!=0) ~ 1 + ratetable(age=age, sex=sex, year=per),
                   ratetable = popm, data = x, method = 'pohar-perme', type = "fleming-harrington", fin.date=ex_date)
  x[, sex := 1L]
  
  ## survtab
  fb <- seq(0, 19, 1/24)
  
  x[, lex.dur := lex.dur/365.242199]
  x[, age := age/365.242199]
  x[, per := get.yrs(per, year.length = "approx")]
  
  setnames(pm, c("year", "agegroup"), c("per", "age"))
  st.e2 <- survtab(Surv(fot, event = lex.Xst) ~ 1, data = x, surv.type="surv.rel", 
                       relsurv.method="e2", pophaz = pm, breaks = list(fot = fb))
  st.pp <- survtab(Surv(fot, event = lex.Xst) ~ 1, data = x, surv.type="surv.rel", 
                       relsurv.method="pp", pophaz = pm, breaks = list(fot = fb))
  setDT(st.e2)
  setDT(st.pp)
  
  ## rs.surv
  fb <- fb[-1]
  fbd <- fb*365.242199
  
  su.e2 <- summary(rs.e2, times = fbd)
  su.e2 <- cbind(data.table(time = fb), data.table(su.e2$surv))
  su.pp <- summary(rs.pp, times = fbd)
  su.pp <- cbind(data.table(time = fb), data.table(su.pp$surv))
  
  expect_equal(st.e2[, r.e2] ,  su.e2[, V1], tolerance = 0.000226, scale = 1L)
  expect_equal(st.pp[, r.pp] ,  su.pp[, V1], tolerance = 0.00292, scale = 1L)
})

# relpois vs. relsurv::rsadd ---------------------------------------------

test_that("relpois congruent with relsurv::rsadd", {
  skip_on_cran()
  
  library(survival)
  library(relsurv)
  
  # male
  pm <- data.table(popEpi::popmort)
  # pm[, surv := 1L]
  pm[, surv := exp(-haz)]
  pm.m <- cast_simple(pm[sex==0], columns = 'year', rows = 'agegroup',  values='surv')
  pm.m[,agegroup := NULL]
  pm.m <- as.matrix(pm.m)
  # female
  pm.f <- cast_simple(pm[sex==1], columns = 'year', rows = 'agegroup',  values='surv')
  pm.f[,agegroup := NULL]
  pm.f <- as.matrix(pm.f)
  
  popm <- transrate(pm.m, pm.f, yearlim = c(1951, 2013), int.length = 1)
  
  pm[, surv := NULL]
  
  sire2 <- copy(sire)
  sire2[, Tstop  := as.integer(ex_date - dg_date)]
  sire2[, dg_age := as.integer(dg_date - bi_date)]
  
  
  sire2[, agegr := cut(dg_age/365.25, breaks = c(0,45,70,Inf))]
  x <- lexpand(sire2, birth  = bi_date, entry = dg_date, exit = ex_date,
               status = status %in% 1:2,
               breaks=list(fot=0:5), pophaz=pm)
  rp <- relpois(x, formula = lex.Xst%in%1:2 ~ -1+FOT+agegr)
  setDT(x)
  setattr(x, "class", c("Lexis", "data.table", "data.frame"))
  sire2[, sex := 2L]
  rs <- relsurv::rsadd(Surv(Tstop, event = status %in% 1:2) ~ ratetable(age=dg_age, sex=sex, year=dg_date) + agegr,
                       int = 0:5,
                       ratetable = popm, data = sire2, method = "glm.poi")
  
  
  expect_equal(coef(rp)[1:5] ,  coef(rs)[3:7], tolerance = 0.055, scale=1, check.attributes=FALSE)
})



test_that("Ederer I expected survival curve agrees with survival::survexp", {
  
  library(survival)
  library(relsurv)
  library(Epi)
  
  # male
  pm <- data.table(popEpi::popmort)
  # pm[, surv := 1L]
  pm[, surv := exp(-haz)]
  pm.m <- cast_simple(pm[sex==0], columns = 'year', rows = 'agegroup',  values='surv')
  pm.m[,agegroup := NULL]
  pm.m <- as.matrix(pm.m)
  # female
  pm.f <- cast_simple(pm[sex==1], columns = 'year', rows = 'agegroup',  values='surv')
  pm.f[,agegroup := NULL]
  pm.f <- as.matrix(pm.f)
  
  popm <- transrate(pm.m, pm.f, yearlim = c(1951, 2013), int.length = 1)
  
  pm[, surv := NULL]
  
  sire2 <- sire[dg_date<ex_date, ]
  set.seed(13)
  sire2 <- sire2[sample(x = 1:.N, size = 500, replace = FALSE)]
  sire2[, Tstop  := as.integer(ex_date - dg_date)]
  sire2[, dg_age := as.integer(dg_date - bi_date)]
  
  x <- Lexis(entry = list(age = dg_age, per = dg_date, fot = 0L),
             exit = list(fot = Tstop), 
             exit.status = as.integer(status %in% 1:2),
             entry.status = 0L, data = sire2)
  setDT(x)
  setattr(x, "class", c("Lexis","data.table", "data.frame"))
  
  ## rs.surv
  ## sex must be coded c(1,2) (male, female)
  x[, sex := 2L]
  
  fb <- seq(0, 19, 1/24)
  su <- survexp(~1, data = x, ratetab = popm, method = "ederer", 
                rmap = list(sex = "female", year = per, age = age),
                times = fb*365.242199)
  
  x[, sex := 1L]
  x[, lex.dur := max(fb)] ## not really needed but illustrative
  x[, age := age/365.242199]
  x[, per := get.yrs(per, year.length = "approx")]
  
  setnames(pm, c("year", "agegroup"), c("per", "age"))
  e1 <- comp_e1(x, breaks = list(fot = fb), pophaz = pm, survScale = "fot")
  
#   plot(su, ylim = c(0.35, 1), col = 1, xscale = 365.242199)
#   lines(surv.exp ~ fot, col = "red", type = "s", data = e1)
  
  ## rs.surv
  fb <- fb[-1]
  fbd <- fb*365.242199
  
  su <- data.table(time = su$time, surv.exp = su$surv)
  su <- su[time != 0L]
  su[, time := time / 365.242199]
  setnames(su, "time", "fot")
  
  expect_equal(e1, su, tolerance = 0.000004575, scale = 1L)
  expect_equal(max(abs(e1$surv.exp-su$surv.exp)), 0L, 
               tolerance = 0.0000103, scale = 1L)
  
  
})

