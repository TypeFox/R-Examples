context("sim_setup")
test_that("sim_setup", {
  tmp <- sim_setup(sim_base(base_id(nDomains = 3, nUnits = 4)) %>% 
                     sim_gen_x() %>%
                     sim_gen_e(), simName = "")
  
  expect_equal(length(tmp), (2))
  expect_that(all(sapply(tmp, inherits, what = "sim_fun")), is_true())
})

test_that("methods equal", {
  setup <- sim_base() %>% sim_gen_x() %>% sim_gen_e() %>% sim_agg()
  cat("\n")
  dat <- show(setup)
  testthat::expect_equal(nrow(dat), 100)
  cat("\n")
  dat <- show(setup %>% sim_comp_agg(function(dat) list(1)))
  testthat::expect_equal(dat, list(1))
  cat("\n")
  dat <- show(setup %>% sim_comp_agg(function(dat) asS4(list(1))))
  testthat::expect_equal(dat, asS4(list(1)))
  cat("\n")
})

context("sim_setup-methods")
test_that("as.data.frame", {
  setup <- sim_base_lm()
  dat <- as.data.frame(setup)
  expect_is(dat, class="data.frame")
  expect_equal(nrow(dat), 10000)
  expect_true(all(names(dat) %in% c("idD", "idU", "y", "x", "e")))
})

test_that("autoplot", {
  expect_warning(autoplot(sim_base_lm(), x = "z"))
  expect_warning(autoplot(sim_base_lm(), y = "k"))
})

test_that("Id construction for not simulated data.frames", {
  dat <- data.frame(id = 1:10, x = rnorm(10))
  dat1 <- sim_base(base_add_id(dat, "id")) %>% as.data.frame
  resultList <- sim_setup(dat, simName = "testthat") %>% sim(R = 10)
  
  expect_equal(length(resultList), 10)
  expect_true(all(names(resultList[[2]]) %in% c("idD", "idU", "idR", "simName", "x", "id")))
})

test_that("sim_setup sorts its content", {
  setup <- sim_base(base_id(3, 3)) %>% sim_agg() %>% sim_sample() %>% sim_comp_n() %>% sim_comp_N() %>% sim_gen_v() %>% sim_gen_x()
  setup1 <- sim_base(base_add_id(data.frame(var = rep(c(1, 2), 10)), "var")) %>% sim_gen_x() %>% sim_gen_e() %>% sim_agg() %>% sim_resp_eq(y = 100 + x + e)
  orderAttr <- sapply(setup, function(fun) fun@order)
  
  expect_equal(length(setup), 6)
  expect_equal(orderAttr, sort(orderAttr))  
  expect_equal(length(setup1), 4)
})