context("CIF's & surv.obs's congruence & comparison w/ survival::survfit")

test_that("surv.obs about the same as Kaplan-Meier & CIFs close to Aalen-Johansen", {
  library(survival)
  BL <- list(fot= seq(0,19,1/12), per=c(2008,2013))
  sire2 <- sire[dg_date<ex_date, ]
  sire2$statusf <- factor(sire2$status, levels = 0:2, 
                         labels = c("alive", "canD", "othD"))
  
  x <- lexpand(sire2, 
               birth  = bi_date, entry = dg_date, exit = ex_date,
               status = statusf,
               breaks=BL)
  st <- survtab(Surv(fot, event = lex.Xst) ~ 1, data = x, surv.type="cif.obs")
  setDT(x)
  setattr(x, "class", c("Lexis", "data.table", "data.frame"))
  setDT(st)
  
  test_that("CIFs and surv.obs sum to 1", {
    expect_equal(st[, CIF_canD + CIF_othD + surv.obs] ,  rep(1, times = st[,.N]), tolerance = 0.0001, scale=1)
  })
  
  x <- lexpand(sire2, 
               birth  = bi_date, entry = dg_date, exit = ex_date,
               status = statusf,
               breaks = BL["per"])
  
  fb <- setdiff(BL$fot, 0)
  su.km  <- survfit(Surv(time=fot, time2=fot+lex.dur, lex.Xst!="alive") ~ 1, data = x)
  su.km  <- summary(su.km, times = fb)
  su.km  <- cbind(data.table(time = su.km$time), data.table(su.km$surv))
  
  su.cif <- survfit(Surv(time=fot, time2=fot+lex.dur, lex.Xst)~1, data=x)
  su.cif <- summary(su.cif, times = fb)
  su.cif <- cbind(data.table(time = su.cif$time), data.table(su.cif$prev))
  
  expect_equal(st[, surv.obs] ,  su.km[, V1], tolerance = 0.0032, scale=1)
  
  expect_equal(su.cif$V1, st$CIF_canD, tolerance = 0.0022, scale=1)
  expect_equal(su.cif$V2, st$CIF_othD, tolerance = 0.0011, scale=1)
})


# custom status var -------------------------------------------------------

test_that("survtab status argument works as expected", {
  skip_on_cran()
  library(survival)
  BL <- list(fot= seq(0,19,1/12), per=c(2008,2013))
  sr <- sire[dg_date < ex_date, ]
  
  sr$statusf <- factor(sr$status, 0:2, labels = c("alive", "canD", "othD"))
  sr$statusb <- as.integer(sr$status %in% 1:2)
  
  st <- NULL
  x <- lexpand(sr, birth  = bi_date, entry = dg_date, 
               exit = ex_date, status = status)
  expect_error({
    suppressWarnings(st <- survtab(Surv(fot, lex.Xst) ~ 1, data = x, surv.type = "surv.obs", 
                      breaks = list(fot = 0:5)))
  }, regexp = "Some status indicators \\(3648 values in total\\) were NA as a result of using Surv\\(\\). Usual suspects: original status variable has NA values, or you have numeric status variable with more than two levels and you did not assign e.g. type = 'mstate' \\(e.g. Surv\\(time = c\\(1,1,1\\), event = c\\(0,1,2\\), type = 'mstate'\\) works\\).")

  st <- NULL
  x <- lexpand(sr, birth  = bi_date, entry = dg_date, 
               exit = ex_date, status = statusf)
  st <- survtab(Surv(fot, lex.Xst) ~ 1, data = x, surv.type = "surv.obs", 
                    breaks = list(fot = 0:5))
  expect_equal(class(st)[1L], "survtab")
  
  st <- NULL
  x <- lexpand(sr, birth  = bi_date, entry = dg_date, 
               exit = ex_date, status = statusb)
  st <- survtab(Surv(fot, lex.Xst) ~ 1, data = x, surv.type = "surv.obs", 
                    breaks = list(fot = 0:5))
  expect_equal(class(st)[1L], "survtab")
  
})


test_that("survtab works with more complicated estimation", {
  library(Epi)
  library(survival)
  
  x <- data.table(popEpi::sire[sire$dg_date < sire$ex_date, ])
  
  ## pretend some are male
  set.seed(1L)
  x$sex <- rbinom(nrow(x), 1, 0.5)
  
  ## period
  x$period <- cut(year(x$dg_date), c(1993,1998,2003,2008,2013), right = F)
  
  # age group
  x$agegr <- cut(x$dg_age, 4)
  
  x$stat <- factor(x$status, levels = 0:2, 
                   labels = c("alive", "canD", "othD"))
  x$enStat <- factor(rep(0L, nrow(x)), levels = 0:2, 
                     labels = c("alive", "canD", "othD"))
  xl <- Lexis(entry = list(FUT = 0, AGE = dg_age, CAL = get.yrs(dg_date)), 
              exit = list(CAL = get.yrs(ex_date)), 
              data = x,
              exit.status = stat, entry.status = enStat, 
              merge = TRUE)
  
  
  ## observed survival
  st1 <- survtab(Surv(time = FUT, event = lex.Xst) ~ factor(sex, 0:1, c("male", "female")) + period 
                     + adjust(agegr), data = xl, 
                     weights = list(agegr = as.numeric(table(x$agegr))),
                     surv.type = "surv.obs",
                     breaks = list(FUT = seq(0, 5, 1/12)))
  
  ag <- splitLexisDT(xl, breaks = seq(0, 5, 1/12), timeScale = "FUT")
  ag[, lex.Cst := as.integer(lex.Cst)]
  ag[, lex.Cst := 0L]
  ag[, lex.Xst := as.integer(lex.Xst != "alive")]
  ag <- aggre(ag, by = list(sex, period, agegr, FUT))
  
  st2 <- survtab_ag(FUT ~ factor(sex, 0:1, c("male", "female")) + period + adjust(agegr), 
                    weights = list(agegr = as.numeric(table(x$agegr))),
                    data = ag, surv.type = "surv.obs", d = "from0to1")

  expect_equal(st1$surv.obs.as.lo, st2$surv.obs.as.lo)
  expect_equivalent(st1, st2)
})



