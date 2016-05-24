
context("Preserve attributes:")
test_that("Attributes are preserved", {
  setup <- base_id(10, 10) %>% sim_gen(addAttr)
  
  expect_equal(attr(as.data.frame(setup), "x"), 1)
  expect_equal(attr(as.data.frame(setup %>% sim_gen_e), "x"), 1)
  expect_equal(attr(as.data.frame(setup %>% sim_gen_v), "x"), 1)
  expect_equal(attr(as.data.frame(setup %>% sim_comp_n), "x"), 1)
  expect_equal(attr(as.data.frame(setup %>% sim_sample(sample_fraction(0.5))), "x"), 1)
  expect_equal(attr(as.data.frame(setup %>% sim_sample(sample_number(2))), "x"), 1)
  expect_equal(attr(as.data.frame(setup %>% sim_sample(sample_numbers(1:10, groupVars = "idD"))), "x"), 1)
  expect_equal(attr(as.data.frame(setup %>% sim_gen_x %>% sim_agg()), "x"), 1)
  expect_equal(attr(as.data.frame(setup %>% sim_resp_eq(y = 1)), "x"), 1)
  expect_equal(attr(as.data.frame(setup %>% sim_gen_ec), "x"), 1)
  
})

