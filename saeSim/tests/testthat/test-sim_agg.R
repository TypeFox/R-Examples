context("sim_agg")

test_that("Basic functionality", {
  setup <- sim_base() %>% sim_gen_e() %>% sim_gen_x() %>% 
    sim_resp_eq(y = 100 + x + e) %>%
    sim_comp_N() %>% 
    sim_comp_pop(comp_var(aFactor = as.factor(c("a", "b") %>% rep(length.out = length(y))),
                          aCharacter = c("a", "b") %>% rep(length.out = length(y)))) %>% 
    sim_agg()
  
  dat <- sim(setup)[[1]]
  
  expect_equal(nrow(dat), 100)
  expect_that(dat$idU, is_a("NULL"))
  expect_that(dat$yMake_factor, is_a("NULL"))
  expect_that(dat$yMake_character, is_a("NULL"))
  expect_that(dat$idD, is_a("numeric"))
  expect_that(dat$e, is_a("numeric"))
  expect_that(dat$x, is_a("numeric"))
  expect_that(dat$y, is_a("numeric"))
  expect_that(dat$N, is_a("numeric"))
  expect_that(dat$aFactora, is_a("numeric"))
  expect_that(dat$aFactorb, is_a("numeric"))
  expect_that(dat$aCharactera, is_a("numeric"))
  expect_that(dat$aCharacterb, is_a("numeric"))
  expect_that(dat$idR, is_a("integer"))
  expect_that(dat$simName, is_a("character"))  
})


test_that("Attributes are preserved", {
  setup <- sim_base() %>% sim_gen_e() %>% sim_gen_x() %>% 
    sim_resp_eq(y = 100 + x + e) %>%
    sim_agg() %>%
    sim_comp_pop(function(dat) {attr(dat, "x") <- 2; dat})
  
  setup %>% as.data.frame %>% attr("x") %>% expect_equal(2)
  
})