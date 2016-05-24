context("Simulator Class")

test_that("simulator_class base class works", {
  sim <- simulator_class$new(100)
  expect_equal(sim$get_name(), "TEMPLATE")
  expect_equal(sim$get_priority(), 100)
  expect_error(sim$simulate())
  expect_error(sim$get_cmd())
  expect_true(is_simulator(sim))
})


test_that("simulator_classs can be registered", {
  register_simulator(simulator_class$new(100))
  expect_true(is_simulator(get_simulator("TEMPLATE")))

  expect_error(register_simulator(1:10))
  expect_error(register_simulator(feature))
})


test_that("reducing of simulation commands works", {
  df <- data.frame(locus_number = rep(1, 10),
                   command = as.character(1:10))
  expect_equal(df, reduce_sim_commands(df))

  df$command <- c(rep("a", 5), rep("b", 5))
  expect_equivalent(reduce_sim_commands(df),
                    data.frame(locus_number = c(5, 5),
                               command = c("a", "b"),
                               stringsAsFactors = FALSE))
})


test_that("simulators can be listed", {
  sims <- list_simulators()
})
