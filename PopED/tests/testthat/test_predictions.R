context("Model Predictions")

test_that("model_prediction works", {
  
  source("examples_fcn_doc/examples_model_prediction.R")
  
  expect_equal(length(unique(df_2$ID)),32)
  
  expect_null(df_3$DV)
  
  expect_null(df_4$a_i)
  
  expect_equal(length(unique(df_5$Group)),2)
  expect_equal(length(unique(df_5$a_i)),2)
  expect_equal(length(unique(df_5$ID)),6)
  
  expect_equal(length(unique(df_6$Group)),2)
  expect_true(all(is.na(df_6$PRED)))
  
  expect_true(all(c("WT","AGE") %in% names(df_7)))
  
  expect_equal(length(unique(df_8$WT)),2)
  expect_equal(length(unique(df_8$AGE)),2)
  
  expect_equal(length(unique(df_9$WT)),2)
  expect_equal(length(unique(df_9$AGE)),2)
  expect_equal(length(unique(df_9$ID)),6)
  
  expect_equal(length(unique(df_10$WT)),6)
  
  expect_equal(length(unique(df_11$AGE)),6)
  
  expect_equal(length(unique(df_12$AMT)),3)
  
  expect_equal(length(unique(df_13$AMT)),2)
  
  expect_equal(length(unique(df_15$AMT[df_15$ID==1])),3)
  
  expect_true("test.csv" %in% list.files())
  
  unlink("test.csv")
  
  dosing_2 <- list(list(AMT=1000,RATE=NA,Time=0.5),list(AMT=3000,RATE=NA,Time=0.5),list(AMT=6000,RATE=NA,Time=0.5))
  
  expect_error(model_prediction(design=design_3,DV=T,dosing=dosing_2))
  
  
})

test_that("plot_model_prediction works", {

  source("examples_fcn_doc/examples_plot_model_prediction.R")
  
})
