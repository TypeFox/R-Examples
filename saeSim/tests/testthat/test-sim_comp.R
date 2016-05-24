context("sim_comp")

test_that("sim_comp and comp_var", {
  setup <- sim_base(base_id(nDomains=5, nUnits = 10)) %>% 
    sim_gen(addAttr) %>%
    sim_gen_x(mean=50, sd=20) %>%
    sim_gen_e(0, 1) %>%
    sim_gen_ec() %>%
    sim_sample(sample_number(size=5, groupVars = "idD")) %>%
    sim_resp_eq(y = 10 * x + e)
  
  dat <- setup %>% 
    sim_comp_popMean() %>% sim_comp_popVar %>% sim_comp_N() %>%
    sim_comp_n() %>% 
    saeSim:::sim_run_once()
  
  calc_lm <- function(dat) {
    linearModel <- lm(y ~ x, data = dat)
    dat$linearPredictor <- predict(linearModel)
    dat
  }
  
  dat1 <- sim(setup %>% sim_agg() %>% sim_comp_agg(calc_lm))[[1]]
  
  expect_equal(class(dat), "data.frame")
  expect_equal(nrow(dat), 25)
  expect_that(all(c("popMean", "popVar", "N", "n") %in% names(dat)), is_true())
  expect_equal(unique(dat$n), 5)
  expect_equal(unique(dat$N), 10)
  expect_that(dat1$linearPredictor, is_equivalent_to(predict(lm(y ~ x, data = dat1))))
  expect_equal(attr(dat, "x"), 1)
  expect_equal(attr(dat1, "x"), 1)
})