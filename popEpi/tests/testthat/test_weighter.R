context("makeWeightsDT")


test_that("makeWeightsDT works with various weights & adjust spesifications", {
  skip_on_cran()
  sibr <- popEpi::sibr[1:100]
  sibr[1:50, sex := 0L]
  
  dt <- lexpand(sibr, birth = bi_date, entry = dg_date, exit = ex_date,
                status = status %in% 1:2, 
                aggre = list(sex, dg_age = popEpi:::cutLow(dg_age, 0:125)))
  
  print <- quote(list(gender = sex))
  prdt <- evalPopArg(dt, print)
  
  adjust <- quote(list(agegr = popEpi:::cutLow(dg_age, c(0, seq(15,85,5), Inf))))
  addt <- evalPopArg(dt, adjust)
  
  expr <- quote(list(pyrs = sum(pyrs), from0to1 = sum(from0to1)))
  dt <- dt[, eval(expr), keyby = cbind(prdt, addt)]
  
  wdt <- ICSS[, list(agegr = c(0, seq(15,85,5)), weights = ICSS3)]
  
  
  expr <- list(quote(list(pyrs = pyrs, from0to1 = from0to1)))
  
  ## works with weights data.frame?
  ## Intention of test: wdt contains more levels of agegr than data;
  ## weights need to be scaled to one correctly and only the original
  ## rows in data kept.
  DT <- makeWeightsDT(dt, values = expr, weights = wdt, 
                      print = "gender", adjust = "agegr")
  expect_equal(nrow(DT), 24L)
  expect_equal(sum(dt$pyrs), sum(DT$pyrs))
  DTW <- unique(DT, by = c("agegr"))$weights
  wdt[, ww := weights/sum(weights)]
  expect_equal(DTW, wdt[agegr %in% DT$agegr, ww])
  expect_equal(DT[, sum(pyrs), keyby = agegr], dt[, sum(pyrs), keyby = agegr], check.attributes = FALSE)
  expect_equal(DT[, sum(pyrs), keyby = gender], dt[, sum(pyrs), keyby = gender], check.attributes = FALSE)
  
  
})



test_that("internal weights are computed correctly", {
  skip_on_cran()
  sibr <- popEpi::sibr[1:100]
  sibr[1:50, sex := 0L]
  
  dt <- lexpand(sibr, birth = bi_date, entry = dg_date, exit = ex_date,
                status = status %in% 1:2, 
                aggre = list(sex, dg_age = popEpi:::cutLow(dg_age, 0:125)))
  
  
  adjust <- quote(popEpi:::cutLow(dg_age, c(0, seq(15,85,5), Inf)))
  dt[, c("agegr") := eval(adjust)]
  
  vals <- list(c("pyrs", "from0to1"))
  
  expect_error({
    DT <- makeWeightsDT(dt, values = vals, 
                        weights = "internal", 
                        print = NULL, adjust = c("sex", "agegr"))
  })
  
  ## one adjust var
  DT1 <- makeWeightsDT(dt, values = vals, 
                       weights = "internal", 
                       print = NULL, adjust = "agegr",
                       internal.weights.values = "pyrs")
  DT2 <- dt[, sum(pyrs)/sum(dt$pyrs), keyby = agegr]  
  expect_equal(DT1$weights, DT2$V1)
  
  ## two adjust vars
  DT1 <- makeWeightsDT(dt, values = vals, 
                       weights = "internal", 
                       print = NULL, adjust = c("sex", "agegr"),
                       internal.weights.values = "pyrs")
  DT2 <- dt[, sum(pyrs), keyby = agegr]  
  DT3 <- dt[, sum(pyrs), keyby = sex] 
  DT4 <- CJ(sex = 0:1, agegr = sort(unique(dt$agegr)))
  DT4 <- merge(DT4, DT3, by = "sex", all.x = TRUE, all.y = TRUE)
  setnames(DT4, "V1", "w.sex")
  DT4 <- merge(DT4, DT2, by = "agegr", all.x = TRUE, all.y = TRUE)
  setnames(DT4, "V1", "w.agegr")
  DT4[, weight := w.agegr*w.sex]
  DT4[, weight := weight/sum(weight)]
  setkeyv(DT4, c("sex", "agegr"))
  setkeyv(DT1, c("sex", "agegr"))
  
  expect_equal(DT1$weights, DT4$weight)
  
})


test_that("weighter works with a list of values arguments", {
  
  skip_on_cran()
  sibr <- popEpi::sibr[1:100]
  sibr[1:50, sex := 0L]
  
  dt <- lexpand(sibr, birth = bi_date, entry = dg_date, exit = ex_date,
                status = status %in% 1:2, 
                aggre = list(sex, dg_age = popEpi:::cutLow(dg_age, 0:125)))
  
  vals <- list(c("pyrs", "from0to1"), quote(from0to0))
  
  DT1 <- makeWeightsDT(dt, values = vals, 
                       print = "sex")
  expect_equivalent(DT1, dt[, lapply(.SD, sum), by = sex, .SDcols = c("pyrs", "from0to1", "from0to0")])
  
  dt$agegr <- cut(dt$dg_age, c(0, 45, 75, Inf))
  DT2 <- makeWeightsDT(dt, values = vals, 
                       print = c("sex", "agegr"))
  expect_equivalent(DT2, dt[, lapply(.SD, sum), by = list(sex,agegr), .SDcols = c("pyrs", "from0to1", "from0to0")])
  
})


test_that("weighted NA checks work", {
  skip_on_cran()
  sibr <- popEpi::sibr[1:100]
  sibr[1:50, sex := 0L]
  
  dt <- lexpand(sibr, birth = bi_date, entry = dg_date, exit = ex_date,
                status = status %in% 1:2, 
                aggre = list(sex, dg_age = popEpi:::cutLow(dg_age, 0:125)))
  dt[, agegr := cut(dg_age, c(0, 45, 75), right = FALSE)]
  vals <- list(c("pyrs", "from0to1"), quote(from0to0))
  
  naTxt <- "A warning message with counts: %%NA_COUNT%%"
  expect_warning({
    DT1 <- makeWeightsDT(dt, values = list(c("pyrs", "from0to1"), quote(from0to0)), 
                         adjust = "agegr", weights = "internal", 
                         internal.weights.values = "pyrs", NA.text = naTxt)
  }, regexp = "A warning message with counts: 15")

})
