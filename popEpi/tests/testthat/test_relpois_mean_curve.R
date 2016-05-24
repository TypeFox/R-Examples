context("test rpcurve vs. survtab congruence")

test_that("rpcurve and survtab e2 are approximately congruent", {
  skip_on_cran()
  
  sire2 <- copy(sire)[dg_date < ex_date, ]
  sire2[, agegr := cut(dg_age, breaks = c(0,45,70,Inf))]
  
  fb <- c(0,3/12,6/12,1:8,10)
  x <- lexpand(sire2, birth  = bi_date, entry = dg_date, exit = ex_date,
               status = status %in% 1:2,
               fot=fb,  pophaz=data.table(popEpi::popmort))
  setDT(x)
  setattr(x, "class", c("Lexis", "data.table", "data.frame"))
  rp <- relpois(x, formula = lex.Xst %in% 1:2 ~ -1 + FOT+agegr)
  mc <- rpcurve(rp)
  
  x$pop.haz <- NULL
  pm <- data.table(popEpi::popmort)
  setnames(pm, c("year", "agegroup"), c("per", "age"))
  w <- as.numeric(table(x$agegr))
  st <- survtab(Surv(fot, lex.Xst) ~ adjust(agegr), 
                pophaz = pm, weights = w,
                relsurv.method = "e2",
                data= x, breaks = list(fot = seq(0, 10, 1/12)))
  setDT(mc)
  setDT(st)
  
  expect_equal(st[Tstop %in% fb]$r.e2.as, mc$est, tolerance=0.0136, scale=1L)
  
  ## added old results on 2016-03-19, 
  ## ref = 4feb1ca37489737332cebf33d24550a9951a7630
  old_res <- c(0.9253749, 0.8801775, 0.8123319, 0.7237591, 0.6679470, 
               0.6315218, 0.6035761, 0.5826368, 0.5645760, 0.5519045, 0.5368186)
  
  expect_equal(old_res, mc$est, tolerance=1e-5, scale=1L)
})



# comparison with flexsurv; maybe not needed ----------------------------------
# library(flexsurv)
# sire2 <- lexpand(sire2, fot=c(0, 10), status)
# sire2[, year := year(ex_date)]
# sire2[, agegroup := as.integer(as.integer(ex_date-bi_date)/365.25)]
# sire2[agegroup > 100, agegroup := 100L]
# 
# sire2 <- data.table:::merge.data.table(sire2, popmort, all.x=F, all.y=F, by=c("sex","year","agegroup"))
# sire2[lex.Xst == 0, haz := 0] ## not really needed, nothing changes even when it works
# 
# ## spline model does not work with bhazard (nothing changes)
# fl <- flexsurvspline(Surv(lex.dur, lex.Xst %in% 1:2) ~ agegr, data=sire2, k = 2, bhazard=sire2$haz)
# ## this works
# fl <- flexsurvreg(Surv(lex.dur, lex.Xst %in% 1:2) ~ agegr, data=sire2, dist="gengamma", bhazard=sire2$haz)
# su <- summary.flexsurvreg(fl, newdata = sire2, ci = FALSE, t = fb, B = 0)
# su <- rbindlist(su)
# su <- su[, list(netsurv = mean(est)), by=time]
# 
# plot(netsurv~time, data=su, type="l")
# lines(est~Tstop, data=mc, col="red")
# # lines(lo~Tstop, data=mc, col="red")
# # lines(hi~Tstop, data=mc, col="red")
# lines(st, "r.e2.as", col="blue", conf.int=FALSE)





