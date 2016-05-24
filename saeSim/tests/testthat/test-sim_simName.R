context("sim_simName")
test_that("sim_simName adds name and returns sim_setup", {
  setup1 <- sim_base_lmmc() 
  setup2 <- setup1 %>% sim_simName("test")
  expect_equal(setup1@simName, "")
  expect_equal(setup2@simName, "test")
  expect_is(setup1, "sim_setup")
})