context("sim_base")
test_that("base_add_id", {
  data(diamonds, envir=environment(), package = "ggplot2")
  diamonds <- diamonds[1:1000, ]
  diamonds["clusterVariable"] <- 1:nrow(diamonds)
  setup <- sim_base(base_add_id(data = diamonds, "clusterVariable"))
  dat <- sim(setup %>% sim_gen_e())[[1]]
  
  expect_equal(nrow(dat), 1000)
  expect_equal(names(dat), c(c("idD", "idU"), names(diamonds), c("e", "idR", "simName")))
})

test_that("sim_base return setup with simName", {
  setup <- sim_base()
  expect_equal(length(setup@simName), 1)
})