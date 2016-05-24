context("%>%")
test_that("basic functionality", {
  setup <- sim_base(base_id(nDomains = 10, nUnits = 100)) %>% sim_gen_e()
  setup1 <- sim_base(base_id(nDomains = 10, nUnits = 10)) %>% sim_gen_x()
  
  expect_equal(dim(sim(setup %>% sim_gen_v())[[1]]), c(1000, 6))
  expect_equal(dim(sim(setup1 %>% sim_gen_ec())[[1]]), c(100, 7))
})
