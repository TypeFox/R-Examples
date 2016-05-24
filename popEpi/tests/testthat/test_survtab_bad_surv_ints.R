context("Testing empty survival intervals in survtab")

test_that("removing consecutively bad surv.ints is logical w/ & w/out adjusting", {
  
  sire2 <- sire[dg_date < ex_date, ]
  sire2[, agegr := cut(dg_age, c(0,45,60,Inf), right=FALSE, labels=FALSE)]
  
  BL <- list(fot= seq(0,10,1/12), per=c(2008,2013))
  
  x <- lexpand(sire2, birth  = bi_date, entry = dg_date, exit = ex_date,
               status = status %in% 1:2, breaks=BL)
  setDT(x)
  
  setattr(x, "class", c("Lexis", "data.table", "data.frame"))
  ## NOTE: neither should give any messages!
  expect_message({
    st1 <-  survtab(Surv(fot, lex.Xst) ~ agegr,
                        subset = !(agegr == 1L & fot > 8.49),
                        data = x, surv.type="surv.obs")
    }, regexp = NA)
  
  ## INTENTION: 7.5+ intervals empty for one age group.
  ## this should make adjusted estimates missing altogether for 7.5+.
  expect_message({
    st2 <- survtab(Surv(fot, lex.Xst) ~ adjust(agegr), 
                      data = x, surv.type="surv.obs",
                      subset = !(agegr == 1L & fot > 8.49),
                      weights = list(agegr = c(0.33, 0.33, 0.33)))
  }, regexp = NA)
  setDT(st1)
  setDT(st2)
  
  
  expect_equal(st1[agegr==3 & Tstop>8.5, .N] ,  18L)
  expect_equal(st1[agegr==1 & Tstop>8.5, .N] ,  0L)
  expect_equal(st2[Tstop > 8.5, .N] , 0L)
})

## non-consecutively bad surv.ints ---------------------------------------------

test_that("survtab_ag messages & results due to non-consecutively bad surv.ints are OK", {
  ## non-consecutively bad surv.ints (missing years 5-6)
  sire2 <- sire[dg_date < ex_date, ]
  sire2[, agegr := cut(dg_age, c(0,45,60,Inf), right=FALSE, labels=FALSE)]
  sire2 <- sire2[!(dg_age > 60 & as.integer(as.integer(ex_date-dg_date)/365.25) %in% 5:6)]
  BL <- list(fot= seq(0,10,1/12), per=c(2008,2013))
  x <- lexpand(sire2, birth  = bi_date, entry = dg_date, exit = ex_date,
               status = status %in% 1:2,
               breaks=BL)
  tf1 <- quote(
    st1 <-  survtab(Surv(fot, lex.Xst)~1, data = x, surv.type="surv.obs", 
                        subset=!(fot >= 5 & fot < 7))
  )
  
  tf2 <- quote(
    st2 <- survtab(Surv(fot, lex.Xst)~adjust(agegr), data = x, surv.type="surv.obs",
                       subset=!(agegr==3 & fot >= 5 & fot < 7), 
                       weights = list(agegr = c(0.33, 0.33, 0.33)))
  )
  
  ## NOTE: \\ needed before "(" or ")"
  msgs <- c(paste0("The total person-time was zero in some survival ",
                   "intervals summed to the margins \\(over any stratifying ",
                   "/ adjusting variables\\) _non-consecutively_, i.e. some ",
                   "intervals after an empty interval had person-time in ",
                   "them. Keeping all survival intervals with some estimates ",
                   "as NA for inspection."),
            "Some cumulative surv.obs were zero or NA:")
  expect_message(eval(tf1), msgs[1],ignore.case=TRUE)
  expect_message(eval(tf1), msgs[2],ignore.case=TRUE)
  
  setDT(st1)
  
  expect_equal(st1[is.na(surv.obs), .N], 60L)
  
  msgs <- c(paste0("The total person-time was zero in some survival ",
                   "intervals, when summed to the variable\\(s\\) ",
                   "'agegr' \\(i.e. over all other variables, if any",
                   "\\) _non-consecutively_, i.e. some intervals after ",
                   "an empty interval had person-time in them. ",
                   "Keeping all survival intervals with some ",
                   "estimates as NA for inspection."),
            "Some cumulative surv.obs were zero or NA:")
  
  expect_message(eval(tf2), msgs[1])
  expect_message(eval(tf2), msgs[2])
  
  setDT(st2)
  expect_equal(st2[is.na(surv.obs.as), .N], 60L)
})





