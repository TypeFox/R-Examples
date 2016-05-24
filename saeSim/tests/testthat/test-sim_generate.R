context("sim_generate")
test_that("sim_generate for smstp_fe", code={
  test_out <- sim_base(base_id(nDomains = 2, nUnits = c(3, 5))) %>%
    sim_gen_x() %>% as.data.frame
  
  expect_is(test_out, "data.frame")
  expect_equal(length(test_out), 3)
  expect_equal(nrow(test_out), 8)
  expect_equal(max(test_out$idU), 5)
  expect_equal(max(test_out$idD), 2)
})

test_that("sim_gen", code={
  setup1 <- sim_base() %>% 
    sim_gen(gen_norm(0, 4, name = "x")) %>% 
    sim_gen(gen_norm(0, 4, "e")) %>% 
    sim_gen_cont(gen_norm(0, 150, "e"), nCont = 0.05, type = "unit", areaVar = "idD", fixed = TRUE)
  setup2 <- sim_base() %>% sim_gen_x() %>% sim_gen_e() %>% sim_gen_ec()
  
  set.seed(1)
  result1 <- sim(setup1, R = 1)
  set.seed(1)
  result2 <- sim(setup2, R = 1)
  
  expect_equal(result2, result1)
})

test_that("gen_generic", {
  # Generator itself
  set.seed(1)
  dat1 <- gen_generic(runif, groupVars = "idD", name = "x")(base_id(5, 2))
  set.seed(1)
  dat2 <- gen_generic(runif, name = "x")(base_id(5, 2))
  expect_is(dat1, "data.frame")
  expect_equal(nrow(dat1), 10)
  expect_equal(dat1[1, "x"], dat1[2, "x"])
  expect_equal(unique(dat2["x"]), dat2["x"])
  expect_error(gen_generic(runif, level = "")(5, 2, "x"))
  
  # In a set-up
  set.seed(1)
  dat1 <- sim(sim_base() %>% 
      sim_gen(gen_generic(rnorm, mean = 0, sd = 4, name="e")) %>%
      sim_gen(gen_generic(rnorm, mean = 0, sd = 1, groupVars = "idD", name="v")))
  set.seed(1)
  dat2 <- sim(sim_base() %>% sim_gen_e() %>% sim_gen_v())
  expect_equal(dat1, dat2)
})

test_that("sim_gen_generic", {
  set.seed(1)
  dat1 <- base_id(5, 5) %>% 
    sim_gen_generic(rnorm, mean = 0, sd = 4, name="e") %>%
    sim_gen_generic(rnorm, mean = 0, sd = 1, groupVars = "idD", name="v") %>% 
    as.data.frame
  
  set.seed(1)
  dat2 <- base_id(5, 5) %>% sim_gen_e() %>% sim_gen_v() %>% as.data.frame
  expect_equal(dat1, dat2)
})

test_that("gen_v_ar1", {
  set.seed(1)
  dat <- base_id_temporal(3, 1, 3) %>%
    sim_gen(gen_v_ar1(
      1.2, sd = 5, rho = 0.6, 
      groupVar = "idD", timeVar = "idT", name = "v")) %>%
    as.data.frame
  testthat::expect_equal(length(unique(dat$v)), 9)
  
  dat <- base_id_temporal(3, 1:3, 3) %>%
    sim_gen(gen_v_ar1(
      1.2, sd = 5, rho = 0.6, 
      groupVar = "idD", timeVar = "idT", name = "v")) %>%
    as.data.frame
  testthat::expect_equal(length(unique(dat$v)), 9)
  
  dat <- base_id_temporal(3, 1:3, 3) %>%
    sim_gen(gen_v_ar1(
      1.2, sd = 5, rho = 0.6, 
      groupVar = c("idD", "idU"), timeVar = "idT", name = "v")) %>%
    as.data.frame
  testthat::expect_equal(length(unique(dat$v)), 18)
  
  # To check that this is somewhat correct...
  # dat <- base_id_temporal(1, 1, 50) %>%
  #   sim_gen(gen_v_ar1(
  #     0, sd = 1, rho = 0.6, 
  #     groupVar = c("idD"), timeVar = "idT", name = "v")) %>%
  #   as.data.frame
  # acf(dat$v)
  
})
