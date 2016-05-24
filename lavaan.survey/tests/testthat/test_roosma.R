

library(lavaan.survey)

context("Roosma example (ess4.gb)")
data(ess4.gb)
  

# Fit the model using lavaan
fit.cfa.ml <- lavaan("range =~ gvjbevn + gvhlthc + gvslvol + gvslvue + gvcldcr + gvpdlwk
   goals =~ sbprvpv  +  sbeqsoc  +  sbcwkfm", data = ess4.gb, estimator = "MLM",
  meanstructure = TRUE, int.ov.free = TRUE, auto.var = TRUE, 
  auto.fix.first = TRUE, auto.cov.lv.x = TRUE)
fit.cfa.ml

# Define the complex survey design for ESS 4 in the UK
des.gb <- svydesign(ids = ~psu, strata = ~stratval, weights = ~dweight, data = ess4.gb)

# Fit the two-factor model while taking the survey design into account.
fit.cfa.surv <- lavaan.survey(fit.cfa.ml, survey.design = des.gb)

fit.cfa.surv.wls.yb <- lavaan.survey(fit.cfa.ml, survey.design=des.gb,
                                      estimator="WLS",
                                      estimator.gamma="Yuan-Bentler")



test_that("scaled chi square matches (MLM)", {
  fm <- fitMeasures(fit.cfa.surv)
  
  expect_equal(fm['chisq.scaling.factor'], 1.540, tolerance = 1e-3, check.attributes = FALSE)
  expect_equal(fm['chisq'], 513.094, tolerance = 1e-3, check.attributes = FALSE)
  expect_equal(fm['df'], 26, check.attributes = FALSE)
  
  expect_equal(fitMeasures(fit.cfa.surv.wls.yb)['chisq'], 269.582, tolerance = 1e-4, check.attributes = FALSE)
})


test_that("standard errors match", {
  skip_on_cran()
  
  ses_surv <- sqrt(diag(vcov(fit.cfa.surv)))
  ses_test_surv <- structure(c(0.0501551740166388, 0.0496972854140471, 0.0467162760651782, 0.0499755056298606, 0.059122658401302, 0.142737184899616, 0.102093768233144, 0.202704445283621, 0.148968353480895, 0.0751907156882466, 0.15348343758362, 0.144747724008307, 0.155749827677067, 0.0300242679984785, 0.0464447623435491, 0.0252751990938959, 0.169859567189274, 0.0299789664048944, 0.0328907149580591, 0.0720970482808726, 0.0407135604060781, 0.0390479224843098, 0.0640747345541657, 0.0585145971200643, 0.0567631513600541, 0.0211856296508052, 0.0258723535785482, 0.0220557069768292), .Names = c("range=~gvhlthc", "range=~gvslvol", "range=~gvslvue", "range=~gvcldcr", "range=~gvpdlwk", "goals=~sbeqsoc", "goals=~sbcwkfm", "gvjbevn~~gvjbevn", "gvhlthc~~gvhlthc", "gvslvol~~gvslvol", "gvslvue~~gvslvue", "gvcldcr~~gvcldcr", "gvpdlwk~~gvpdlwk", "sbprvpv~~sbprvpv", "sbeqsoc~~sbeqsoc", "sbcwkfm~~sbcwkfm", "range~~range", "goals~~goals", "range~~goals", "gvjbevn~1", "gvhlthc~1", "gvslvol~1", "gvslvue~1", "gvcldcr~1", "gvpdlwk~1", "sbprvpv~1", "sbeqsoc~1", "sbcwkfm~1"))

  expect_equal(ses_surv, ses_test_surv, tolerance = 1e-6)
  
  ses_wls <- sqrt(diag(vcov(fit.cfa.surv.wls.yb)))
  ses_test_wls <- structure(c(0.0382872099011375, 0.0380051964453317, 0.0436503524754807, 0.045606888959657, 0.0500624097356357, 0.161870799318694, 0.0916039502138099, 0.211160937898446, 0.117597219961549, 0.0617223783648685, 0.133705333782911, 0.112238961243799, 0.132429614454983, 0.0267461440579609, 0.0439271663666294, 0.0202694631968092, 0.202959121184718, 0.0271010276867731, 0.0349421771706162, 0.0720970482808728, 0.040713560406078, 0.0390479224843097, 0.0640747345541659, 0.0585145971200642, 0.0567631513600539, 0.0211856296508052, 0.0258723535785482, 0.0220557069768292), .Names = c("range=~gvhlthc", "range=~gvslvol", "range=~gvslvue", "range=~gvcldcr", "range=~gvpdlwk", "goals=~sbeqsoc", "goals=~sbcwkfm", "gvjbevn~~gvjbevn", "gvhlthc~~gvhlthc", "gvslvol~~gvslvol", "gvslvue~~gvslvue", "gvcldcr~~gvcldcr", "gvpdlwk~~gvpdlwk", "sbprvpv~~sbprvpv", "sbeqsoc~~sbeqsoc", "sbcwkfm~~sbcwkfm", "range~~range", "goals~~goals", "range~~goals", "gvjbevn~1", "gvhlthc~1", "gvslvol~1", "gvslvue~1", "gvcldcr~1", "gvpdlwk~1", "sbprvpv~1", "sbeqsoc~1", "sbcwkfm~1"))
    
  expect_equal(ses_wls, ses_test_wls, tolerance = 1e-6)

})

test_that("an estimate matches", {
  
  expect_equal(coef(fit.cfa.surv)['range~~range'], 1.892812, tolerance = 1e-5, check.attributes = FALSE)
  expect_equal(coef(fit.cfa.surv.wls.yb)['range~~range'], 2.369951, tolerance = 1e-5, check.attributes = FALSE) 
})

test_that("rep weights works", {
  skip_on_cran()
  
  des.gb.rep <- as.svrepdesign(des.gb, type="JKn")
  fit.cfa.surv.rep <- lavaan.survey(fit.cfa.ml, survey.design=des.gb.rep)
  expect_equal(fitMeasures(fit.cfa.surv.rep)['chisq.scaled'],  332.2483, tolerance = 1e-4, check.attributes = FALSE)

  expect_equal(coef(fit.cfa.surv.rep)['range~~range'], 1.892877, tolerance = 1e-5, check.attributes = FALSE)
  
  ses <- sqrt(diag(vcov(fit.cfa.surv.rep)))
  ses_exp <- structure(c(0.0502666632522885, 0.0498131899551345, 0.0467907385720135, 0.0500510136357558, 0.0592838936370152, 0.142966126809878, 0.102231779136298, 0.202889310870968, 0.149222825840217, 0.0752050895793933, 0.153771077584903, 0.144917328162438, 0.155883264607703, 0.0300566754283666, 0.0465486121901415, 0.025313323826504, 0.170458795096919, 0.0300465605774563, 0.0330329932040509, 0.0721448001181756, 0.0407449049500816, 0.0390678099978025, 0.0641233249663045, 0.058546582138815, 0.0568040171874358, 0.0211966461510079, 0.0258955069765785, 0.0220749866699128), .Names = c("range=~gvhlthc", "range=~gvslvol", "range=~gvslvue", "range=~gvcldcr", "range=~gvpdlwk", "goals=~sbeqsoc", "goals=~sbcwkfm", "gvjbevn~~gvjbevn", "gvhlthc~~gvhlthc", "gvslvol~~gvslvol", "gvslvue~~gvslvue", "gvcldcr~~gvcldcr", "gvpdlwk~~gvpdlwk", "sbprvpv~~sbprvpv", "sbeqsoc~~sbeqsoc", "sbcwkfm~~sbcwkfm", "range~~range", "goals~~goals", "range~~goals", "gvjbevn~1", "gvhlthc~1", "gvslvol~1", "gvslvue~1", "gvcldcr~1", "gvpdlwk~1", "sbprvpv~1", "sbeqsoc~1", "sbcwkfm~1"))
  
  expect_equal(ses, ses_exp)
})