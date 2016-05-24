context("sim_sample")
test_that("Attributes are preserved", {
  setup <- sim_base(base_id(nDomains=3, nUnits = 4)) %>%
    sim_gen_x(mean=50, sd=20, name = "x") %>%
    sim_gen_e(0, 1, name = "e") %>% 
    sim_comp_pop(function(dat) {attr(dat, "x") <- 2; dat})
  
  setup %>% sim_sample(sample_number(2)) %>% as.data.frame %>% attr("x") %>% expect_equal(2)
  setup %>% sim_sample(sample_fraction(0.2)) %>% as.data.frame %>% attr("x") %>% expect_equal(2)
  
  setup %>% sim_sample(sample_number(2, groupVars = "idD")) %>% as.data.frame %>% attr("x") %>% expect_equal(2)
  setup %>% sim_sample(sample_fraction(0.2, groupVars = "idD")) %>% as.data.frame %>% attr("x") %>% expect_equal(2)
  
})

test_that("Basic sampling functionality", {
  setup <- sim_base(base_id(nDomains=3, nUnits = 4))
  
  setup %>% sim_sample(sample_fraction(size = 0.05)) %>% as.data.frame %>% 
    nrow %>% expect_equal(1)
  
  setup %>% sim_sample(sample_number(size = 5L)) %>% as.data.frame %>% 
    nrow %>% expect_equal(5)
  
  setup %>% sim_sample(sample_number(size = 2L, groupVars = "idD")) %>% as.data.frame %>% 
    nrow %>% expect_equal(6)
  
  setup <- base_id(3, 100)
  setup %>% sim_sample(sample_fraction(size = 0.01, groupVars = "idD")) %>% as.data.frame %>% 
    nrow %>% expect_equal(3)
})

test_that("applying the sampling functions correctly", {
  setup <- sim_base(base_id(nDomains=3, nUnits = 10)) %>%
    sim_gen_x(mean=50, sd=20, name = "x")
  
  expect_equal(nrow(sim(setup %>% sim_sample())[[1]]), 15)
  expect_equal(nrow(sim(setup %>% sim_sample(sample_fraction(0.05)))[[1]]), (2))
  
  expect_equal(nrow(sim(setup %>% sim_sample(sample_number(5L)))[[1]]), (5))
  expect_equal(nrow(sim(setup %>% sim_sample(sample_number(c(3), groupVars = "idD")))[[1]]), (9))
})

test_that("more sampling", {
  setup <- sim_base()
  result1 <- setup %>% sim_sample(sample_number(1)) %>% as.data.frame
  result2 <- setup %>% sim_sample(sample_fraction(0.5)) %>% as.data.frame
  result3 <- setup %>% sim_sample(sample_number(5.5)) %>% as.data.frame
  result4 <- setup %>% sim_sample(sample_number(1L)) %>% as.data.frame
  expect_equal(nrow(result1), (1))
  expect_equal(nrow(result2), (5000))
  expect_equal(nrow(result3), (5))
  expect_equal(nrow(result4), (1))
})

test_that("vectorized versions of sample_number", {
  setup <- sim_base()
  result1 <- setup %>% sim_sample(sample_numbers(1, groupVars = "idD")) %>% as.data.frame
  result2 <- setup %>% sim_sample(sample_numbers(rep(c(1, 2), 50), groupVars = "idD")) %>% as.data.frame  
  expect_equal(nrow(result1), 100)
  expect_equal(nrow(result2), 150)
})

test_that("sampling of clusters", {
  dat <- addAttr(base_id())
  
  dat1 <- sample_cluster_number(2, groupVars = "idD")(dat)
  expect_is(dat1, "data.frame")
  expect_equal(attr(dat1, "x"), 1)
  expect_equal(nrow(dat1), 20)
  expect_equal(length(unique(dat1$idD)), 2)
  
  dat2 <- sample_cluster_fraction(0.2, groupVars = "idD")(dat)
  expect_is(dat2, "data.frame")
  expect_equal(attr(dat2, "x"), 1)
  expect_equal(nrow(dat2), 20)
  expect_equal(length(unique(dat2$idD)), 2)
})