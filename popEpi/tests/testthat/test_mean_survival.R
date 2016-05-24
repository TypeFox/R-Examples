context("mean survival testing")

test_that("survmean() agrees with old results", {
  skip_on_cran()
  library(Epi)
  library(survival)
  
  sr <- copy(sire)[dg_date < ex_date, ]
  sr$agegr <- cut(sr$dg_age, c(0,45,60,Inf), right=FALSE)
  
  x <- Lexis(entry = list(FUT = 0, AGE = dg_age, CAL = get.yrs(dg_date)),
             exit = list(CAL = get.yrs(ex_date)),
             data = sr,
             exit.status = factor(status, levels = 0:2,
                                  labels = c("alive", "canD", "othD")),
             entry.status = factor(0, levels = 0:2,
                                   labels = c("alive", "canD", "othD")),
             merge = TRUE)
  
  ## observed survival
  pm <- copy(popEpi::popmort)
  names(pm) <- c("sex", "CAL", "AGE", "haz")
  sm <- survmean(Surv(time = FUT, event = lex.Xst != "alive") ~ agegr,
                 pophaz = pm, data = x,
                 breaks = list(FUT = seq(0, 10, 1/12)),
                 e1.breaks = list(FUT = c(seq(0, 10, 1/12), 11:100)))
  
  ## values to test against computed on 2016-03-04;
  ## git ref: 5077677
  expect_equal(sm$est, c(33.951439, 21.611419,  7.604318), tol = 0.005, scale = 1)
  expect_equal(sm$exp, c(45.25686, 31.22712, 13.06725), tol = 0.005, scale = 1)
  
  
})



test_that("survmean() agrees with results computed using pkg survival", {
  skip_on_cran()
  ## will only compute the mean survival time based on the cohort.
  ## will also compute expected extrapolation curve by hand.
  as.date.Date <- function(x, ...) {
    x <- as.integer(x) + 3653L
    as.date(x)
  }
  #### compute observed survivals
  library(survival)
  library(Epi)
  library(relsurv)
  
  BL <- list(fot= seq(0,15,1/12))
  eBL <- list(fot = unique(c(BL$fot, seq(15, 115,0.5))))
  sire2 <- sire[dg_date<ex_date, ]
  sire2$statusf <- factor(sire2$status, levels = 0:2, 
                          labels = c("alive", "canD", "othD"))
  
  x <- lexpand(sire2, 
               birth  = bi_date, entry = dg_date, exit = ex_date,
               status = statusf,
               breaks = NULL)
  popmort_sm <- setDT(copy(popEpi::popmort))
  setnames(popmort_sm, c("agegroup", "year"), c("age", "per"))
  sm <- survmean(Surv(fot, event = lex.Xst) ~ 1, 
                 breaks = BL, e1.breaks = eBL,
                 pophaz = popmort_sm, data = x)
  st <- survtab(Surv(fot, event = lex.Xst) ~ 1, 
                    breaks = BL,
                    data = x, surv.type="surv.obs")
  setDT(x)
  setattr(x, "class", c("Lexis", "data.table", "data.frame"))
  setDT(st)
  
  fb <- setdiff(BL$fot, 0)
  su.km  <- survfit(Surv(time=fot, time2=fot+lex.dur, lex.Xst!="alive") ~ 1, data = x)
  su.km  <- summary(su.km, times = fb)
  su.km  <- cbind(data.table(time = su.km$time), data.table(su.km$surv))
  
  #### compute expected survivals for all subjects
  BL <- list(fot = seq(0,100, 1/2))
  xe <- copy(x)
  setattr(xe$age, "class", c("yrs", "numeric"))
  setattr(xe$per, "class", c("yrs", "numeric"))
  setattr(xe$per, "year.length", "approx")
  setattr(xe$age, "year.length", "approx")
  xe[, perdate := as.date.Date(as.Date.yrs(per))]
  xe[, agedate := as.integer(as.Date.yrs(per)-bi_date)]
  
  
  ## form ratetable
  # male
  pm <- setDT(copy(popEpi::popmort))
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
  
  ## survival::survexp()
  su.exp <- survexp(~1, data = xe, ratetab = popm, method = "ederer", 
                    rmap = list(sex = "female", year = perdate, 
                                age = agedate),
                    times = BL$fot*365.242199)
  
  su.exp <- data.table(time = su.exp$time/365.242199, surv.exp = su.exp$surv)
  su.exp <- su.exp[time != 0L]
  su.exp[, time := time + max(su.km$time)]
  
  ## popEpi:::comp_e1()
  setnames(pm, c("year", "agegroup"), c("per", "age"))
  empty_list <- list(fot = NULL, per = NULL, age = NULL)
  forceLexisDT(xe, breaks = empty_list, allScales = c("fot", "per", "age"))
  e1 <- comp_e1(xe, breaks = BL, pophaz = pm, survScale = "fot")
  
  # plot(I(c(1,surv.exp)) ~ I(BL$fot), ylim = c(0, 1), col = 1, data = su.exp, type = "s")
  # lines(I(c(1,surv.exp)) ~ I(BL$fot), col = "red", type = "s", data = e1)
  
  su.exp[, surv.exp := surv.exp * su.km[.N, V1]]
  su <- rbindlist(list(su.km, su.exp))
  
  st <- st[, .(fot = Tstop, surv.obs)]
  st <- rbindlist(list(st, e1))
  st[181:380, fot := fot + st[180, fot]]
  st[181:380, surv.obs := surv.obs*st[180L, surv.obs]]
  
  # plot(V1 ~ time, ylim = c(0, 1), col = 1, data = su, type = "l")
  # lines(surv.obs ~ fot, col = "red", type = "l", data = st)
  
  expect_equal(st$surv.obs, su$V1, scale = 1L, tolerance = 0.0001765)
  expect_equal(max(abs(st$surv.obs-su$V1)), 0.0003814132, scale = 1L)
  
  st[, delta := fot - c(0, fot[-.N])]
  st[, l1 := c(1, surv.obs[-.N])]
  sm.st <- st[, sum((surv.obs+l1)/2L*delta)]
  
  su[, delta := time - c(0, time[-.N])]
  su[, l1 := c(1, V1[-.N])]
  sm.su <- su[, sum((V1+l1)/2L*delta)]
  
  expect_equal(sm.st, sm.su, scale = 1L, tolerance = 0.0571)
  expect_equal(abs(sm.st-sm.su), 0.01086865, scale = 1L)
})


test_that("survmean expected survival curve corresponds to full Ederer I", {
  skip_on_cran()
  library(Epi)
  library(survival)
  
  sr <- copy(sire)[dg_date < ex_date, ]
  sr$agegr <- cut(sr$dg_age, c(0,45,60,Inf), right=FALSE)
  
  x <- Lexis(entry = list(FUT = 0, AGE = dg_age, CAL = get.yrs(dg_date)),
             exit = list(CAL = get.yrs(ex_date)),
             data = sr,
             exit.status = factor(status, levels = 0:2,
                                  labels = c("alive", "canD", "othD")),
             entry.status = factor(0, levels = 0:2,
                                   labels = c("alive", "canD", "othD")),
             merge = TRUE)
  
  pm <- copy(popEpi::popmort)
  names(pm) <- c("sex", "CAL", "AGE", "haz")
  
  BL <- list(FUT = seq(0, 10, 1/12))
  eBL <- list(FUT = c(BL$FUT, seq(11,110,1/2)))
  sm <- survmean(Surv(time = FUT, event = lex.Xst != "alive") ~ 1,
                 pophaz = pm, data = x,
                 breaks = BL, 
                 e1.breaks = eBL)
  
  ## pure Ederer I curve
  setDT(x)
  x[, lex.dur := 110]
  setattr(x, "class", c("Lexis", "data.table", "data.frame"))
  e1 <- comp_e1(x = x, breaks = eBL, 
                pophaz = pm, survScale = "FUT")
  setkeyv(e1, "FUT")
  e1[, delta := diff(eBL$FUT)]
  e1[, l1 := c(1,surv.exp[-.N])]
  sm.e1 <- e1[, sum((l1+surv.exp)/2L*delta)]
  
  
  expect_equal(sm$exp, sm.e1, tol = 0.0005, scale = 1)
  
})

test_that("survmean period method is useful", {
  skip_on_cran()
  library(Epi)
  library(survival)
  
  sr <- copy(sire)[dg_date < ex_date, ]
  sr$agegr <- cut(sr$dg_age, c(0,45,60,Inf), right=FALSE)
  
  x <- Lexis(entry = list(FUT = 0, AGE = dg_age, CAL = get.yrs(dg_date)),
             exit = list(CAL = get.yrs(ex_date)),
             data = sr,
             exit.status = factor(status, levels = 0:2,
                                  labels = c("alive", "canD", "othD")),
             entry.status = factor(0, levels = 0:2,
                                   labels = c("alive", "canD", "othD")),
             merge = TRUE)
  
  pm <- data.table(popEpi::popmort)
  names(pm) <- c("sex", "CAL", "AGE", "haz")
  
  BL <- list(FUT = seq(0, 10, 1/12))
  eBL <- list(FUT = c(BL$FUT, seq(max(BL$FUT),110,1/2)))
  eBL$FUT <- sort(unique(eBL$FUT))
  sm <- survmean(Surv(time = FUT, event = lex.Xst != "alive") ~ 1,
                 subset = dg_date >= "1998-01-01" & dg_date < "2003-01-01",
                 pophaz = pm, data = x,
                 r = 1,
                 breaks = BL, 
                 e1.breaks = eBL)
  
  
  BL <- list(FUT = seq(0, 5, 1/12), CAL = c(1998,2003))
  eBL <- list(FUT = c(BL$FUT, seq(max(BL$FUT), 10, 1/12), seq(10,110,1/2)))
  eBL$FUT <- sort(unique(eBL$FUT))
  smp <- survmean(Surv(time = FUT, event = lex.Xst != "alive") ~ 1,
                  pophaz = pm, data = x,
                  breaks = BL, 
                  r = 1,
                  e1.breaks = eBL)
  
  
  
  expect_equal(sm$obs, smp$obs)
  expect_equal(sm$exp, smp$exp)
  expect_equal(smp$est, 10.01542, tol = 0.0005, scale = 1)
  expect_equal(sm$est, 10.0216, tol = 0.0005, scale = 1)
  
})



test_that("Dates and frac. yrs produce congruent results", {
  skip_on_cran()
  library(Epi)
  library(survival)
  
  x <- data.table(popEpi::sire)
  x <- x[dg_date<ex_date]
  
  ## phony group variable
  set.seed(1L)
  x$group <- rbinom(nrow(x), 1, 0.5)
  
  
  ## yrs
  xy <- Lexis(entry = list(FUT = 0, AGE = dg_age, CAL = get.yrs(dg_date)), 
              exit = list(CAL = get.yrs(ex_date)), 
              data = x,
              exit.status = factor(status, levels = 0:2, 
                                   labels = c("alive", "canD", "othD")), 
              entry.status = factor(0, levels = 0:2,
                                    labels = c("alive", "canD", "othD")),
              merge = TRUE)
  
  ## dates
  xd <- Lexis(entry = list(FUT = 0L, AGE = dg_date-bi_date, CAL = dg_date),
              exit = list(CAL = ex_date),
              data = x,
              exit.status = factor(status, levels = 0:2, 
                                   labels = c("alive", "canD", "othD")), 
              entry.status = factor(0, levels = 0:2,
                                    labels = c("alive", "canD", "othD")),
              merge = TRUE)
  yd <- 365.242199
  BLy <- list(FUT = seq(0, 9, 1/4))
  BLd <- lapply(BLy, function(el) el * yd)
  eBLy <- list(FUT = c(seq(0, 9, 1/4), 9.5, 10:75))
  eBLd <- lapply(eBLy, function(el) el * yd)
  
  pmy <- data.table(popEpi::popmort)
  setnames(pmy, c("year", "agegroup"), c("CAL", "AGE"))
  
  pmd <- data.table(pmy)
  pmd[, CAL := as.Date(paste0(CAL, "-01-01"))]
  pmd[, AGE := AGE * yd]
  pmd[, haz := haz/yd]
  
  #### hazard method
  ## observed survival & Ederer II
  
  sty <- survmean(Surv(FUT, lex.Xst) ~ group, data = xy, 
                  surv.method = "hazard",
                  e1.breaks = eBLy,
                  breaks = BLy, pophaz = pmy)
  
  std <- survmean(Surv(FUT, lex.Xst) ~ group, data = xd, 
                  surv.method = "hazard",
                  e1.breaks = eBLd,
                  breaks = BLd, pophaz = pmd)    
  cuy <- data.table(attributes(sty)$survmean.meta$curves)
  cud <- data.table(attributes(std)$survmean.meta$curves)
  
  std[, c("est", "exp", "YPLL") := lapply(.SD, function(col) col/yd),
      .SDcols = c("est", "exp", "YPLL")]
  
  expect_equal(sty$est, std$est, scale = 1L, tolerance = 0.0005)
  expect_equal(sty$exp, std$exp, scale = 1L, tolerance = 0.001)
  expect_equal(sty$obs, std$obs)
  
  expect_equal(cuy$surv, cud$surv, scale = 1L, tolerance = 0.00005)
  
  #### lifetable method
  ## observed survival & Ederer II
  
  sty <- survmean(Surv(FUT, lex.Xst) ~ group, data = xy, 
                  surv.method = "lifetable",
                  e1.breaks = eBLy,
                  breaks = BLy, pophaz = pmy)
  
  std <- survmean(Surv(FUT, lex.Xst) ~ group, data = xd, 
                  surv.method = "lifetable",
                  e1.breaks = eBLd,
                  breaks = BLd, pophaz = pmd)    
  cuy <- data.table(attributes(sty)$survmean.meta$curves)
  cud <- data.table(attributes(std)$survmean.meta$curves)
  
  std[, c("est", "exp", "YPLL") := lapply(.SD, function(col) col/yd),
      .SDcols = c("est", "exp", "YPLL")]
  
  expect_equal(sty$est, std$est, scale = 1L, tolerance = 0.0005)
  expect_equal(sty$exp, std$exp, scale = 1L, tolerance = 0.001)
  expect_equal(sty$obs, std$obs)
  
  expect_equal(cuy$surv, cud$surv, scale = 1L, tolerance = 0.00005)
  
  
})







