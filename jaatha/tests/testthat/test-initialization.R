context("Initialization")

test_that("determining block per par works", {
  expect_equal(determine_bpp(2, 1), 1)
  expect_equal(determine_bpp(2, 2), 2)
  expect_equal(determine_bpp(2, 3), 2)
  expect_equal(determine_bpp(2, 4), 2)
  expect_equal(determine_bpp(2, 5), 3)
  expect_equal(determine_bpp(2, 9), 5)
  
  expect_equal(determine_bpp(3, 1), 1)
  expect_equal(determine_bpp(3, 2), 2)
  expect_equal(determine_bpp(3, 3), 2)
  expect_equal(determine_bpp(3, 4), 2)
  expect_equal(determine_bpp(3, 5), 2)
  expect_equal(determine_bpp(3, 9), 3)
})


test_that("creation of initial blocks works", {
  par_ranges <- par_ranges_class$new(matrix(1:4, 2, 2))
  expect_equal(1, length(create_initial_blocks(par_ranges, 1)))
  expect_equal(4, length(create_initial_blocks(par_ranges, 2)))
  expect_equal(6, length(create_initial_blocks(par_ranges, 3)))
  expect_equal(8, length(create_initial_blocks(par_ranges, 4)))
  par_ranges <- par_ranges_class$new(matrix(1:6, 3, 2))
  expect_equal(1, length(create_initial_blocks(par_ranges, 1)))
  expect_equal(6, length(create_initial_blocks(par_ranges, 2)))
  expect_equal(9, length(create_initial_blocks(par_ranges, 3)))
  expect_equal(12, length(create_initial_blocks(par_ranges, 4)))
})


test_that("a complete initial search works", {
  model <- create_test_model()
  data <- create_test_data(model)
  
  sim_cache <- create_sim_cache()
  par <- do_initial_search(model, data, 1, sim = 20, cores = 1, sim_cache)
  expect_that(par, is_a("matrix"))
  expect_equal(dim(par), c(1, model$get_par_number()))
  expect_true(all(par >= 0 & par <= 1))
  
  sim_cache <- create_sim_cache()
  par <- do_initial_search(model, data, 3, sim = 20, cores = 1, sim_cache)
  expect_that(par, is_a("matrix"))
  expect_equal(dim(par), c(3, model$get_par_number()))
  expect_true(all(par >= 0 & par <= 1))
  
  sim_cache <- create_sim_cache()
  par <- do_initial_search(model, data, 4, sim = 20, cores = 1, sim_cache)
  expect_that(par, is_a("matrix"))
  expect_equal(dim(par), c(4, model$get_par_number()))
  expect_true(all(par >= 0 & par <= 1))
})


test_that("zoom-in search works", {
  model <- create_test_model()
  data <- create_test_data(model)
  
  sim_cache <- create_sim_cache()
  par <- do_zoom_in_search(model, data, 1, sim = 20, cores = 1, sim_cache, 0.05)
  expect_that(par, is_a("matrix"))
  expect_equal(dim(par), c(1, model$get_par_number()))
  expect_true(all(par >= 0 & par <= 1))
  
  sim_cache <- create_sim_cache()
  par <- do_zoom_in_search(model, data, 2, sim = 20, cores = 1, sim_cache, 0.1)
  expect_that(par, is_a("matrix"))
  expect_equal(dim(par), c(2, model$get_par_number()))
  expect_true(all(par >= 0 & par <= 1))
  
  sim_cache <- create_sim_cache()
  par <- do_zoom_in_search(model, data, 3, sim = 20, cores = 1, sim_cache, 0.2)
  expect_that(par, is_a("matrix"))
  expect_equal(dim(par), c(3, model$get_par_number()))
  expect_true(all(par >= 0 & par <= 1))
})


test_that("getting the start positions works", {
  model <- create_test_model()
  data <- create_test_data(model)
  
  # middle
  sim_cache <- create_sim_cache()
  expect_equal(get_start_pos(model, data, 1, 20, "middle", 1, sim_cache),
               matrix(0.5, 1, model$get_par_number()))
  expect_equal(get_start_pos(model, data, 2, 20, "middle", 1, sim_cache),
               matrix(0.5, 2, model$get_par_number()))
  expect_equal(get_start_pos(model, data, 3, 20, "middle", 1, sim_cache),
               matrix(0.5, 3, model$get_par_number()))
  
  # random
  pos <- get_start_pos(model, data, 1, 20, "random", 1, sim_cache)
  expect_equal(dim(pos), c(1, 2))
  expect_equal(length(unique(pos)), 2)
  expect_true(all(pos >= 0 & pos <= 1))
  
  pos <- get_start_pos(model, data, 2, 20, "random", 1, sim_cache)
  expect_equal(dim(pos), c(2, 2))
  expect_equal(length(unique(pos)), 4)
  expect_true(all(pos >= 0 & pos <= 1))
  
  pos <- get_start_pos(model, data, 3, 20, "random", 1, sim_cache)
  expect_equal(dim(pos), c(3, 2))
  expect_equal(length(unique(pos)), 6)
  expect_true(all(pos >= 0 & pos <= 1))
  
  
  # initial search
  sim_cache <- create_sim_cache()
  pos <- get_start_pos(model, data, 1, 20, "initial-search", 1, sim_cache, 0.05)
  expect_that(pos, is_a("matrix"))
  expect_true(all(pos >= 0 & pos <= 1))
  
  # zoom-in
  sim_cache <- create_sim_cache()
  pos <- get_start_pos(model, data, 1, 20, "zoom-in", 1, sim_cache, 0.05)
  expect_that(pos, is_a("matrix"))
  expect_true(all(pos >= 0 & pos <= 1))
  
  # errors
  sim_cache <- create_sim_cache()
  expect_error(get_start_pos(model, data, 1, 20, "1", 1, sim_cache, 0.05))
  expect_error(get_start_pos(model, data, 1, 20, 1, 1, sim_cache, 0.1))
})
