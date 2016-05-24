context("SumStat Trees")

test_that("tree summary statistics require files", {
  expect_false(requires_files(coal_model(5)))
  expect_false(requires_segsites(coal_model(5)))
  expect_true(requires_trees((coal_model(5) + sumstat_trees())))
  expect_true(requires_trees((coal_model(5) + sumstat_sg_trees())))
  expect_false(requires_segsites((coal_model(5) + sumstat_trees())))
  expect_false(requires_segsites((coal_model(5) + sumstat_sg_trees())))
})


test_that("Generating trees for trios works", {
  trees <- list(c("[2](2:0.865,(1:0.015,3:0.015):0.850);",
                  "[3](2:0.865,(1:0.015,3:0.015):0.850);",
                  "[4](2:1.261,(1:0.015,3:0.015):1.246);",
                  "[11](2:1.261,(1:0.015,3:0.015):1.246);"),
                c("[2](3:0.613,(1:0.076,2:0.076):0.537);",
                  "[18](3:0.460,(1:0.076,2:0.076):0.384);"))

  trio_trees <- generate_trio_trees(trees, matrix(c(2, 4, 8, 2, 4, 2), 1, 6))
  expect_that(trio_trees, is_a("list"))
  expect_equal(length(trio_trees), 1)

  expect_equal(trio_trees[[1]][[1]][1], "[2](2:0.865,(1:0.015,3:0.015):0.850);")
  expect_equal(trio_trees[[1]][[2]][1], "[3](2:1.261,(1:0.015,3:0.015):1.246);")
  expect_equal(trio_trees[[1]][[2]][2], "[5](2:1.261,(1:0.015,3:0.015):1.246);")
  expect_equal(trio_trees[[1]][[3]][1], "[4](2:1.261,(1:0.015,3:0.015):1.246);")

  expect_equal(trio_trees[[1]][[1]][2], "[2](3:0.613,(1:0.076,2:0.076):0.537);")
  expect_equal(trio_trees[[1]][[2]][3], "[8](3:0.460,(1:0.076,2:0.076):0.384);")
  expect_equal(trio_trees[[1]][[3]][2], "[4](3:0.460,(1:0.076,2:0.076):0.384);")


  trio_trees_2 <- generate_trio_trees(list(trees[[1]], trees[[2]],
                                           trees[[1]], trees[[2]]),
                                      matrix(c(2, 4, 8, 2, 4, 2), 2, 6, TRUE))
  expect_equal(trio_trees_2, list(trio_trees[[1]], trio_trees[[1]]))

  trio_trees <- generate_trio_trees(trees, matrix(c(9, 2, 2, 5, 2, 2), 1, 6)) #nolint
  expect_equal(trio_trees[[1]][[1]][1], "[2](2:0.865,(1:0.015,3:0.015):0.850);")
  expect_equal(trio_trees[[1]][[1]][2], "[3](2:0.865,(1:0.015,3:0.015):0.850);")
  expect_equal(trio_trees[[1]][[1]][3], "[4](2:1.261,(1:0.015,3:0.015):1.246);")
  expect_equal(trio_trees[[1]][[2]][1], "[2](2:1.261,(1:0.015,3:0.015):1.246);")
  expect_equal(trio_trees[[1]][[3]][1], "[2](2:1.261,(1:0.015,3:0.015):1.246);")
  expect_equal(trio_trees[[1]][[1]][4], "[2](3:0.613,(1:0.076,2:0.076):0.537);")
  expect_equal(trio_trees[[1]][[1]][5], "[7](3:0.460,(1:0.076,2:0.076):0.384);")
  expect_equal(trio_trees[[1]][[2]][2], "[2](3:0.460,(1:0.076,2:0.076):0.384);")
  expect_equal(trio_trees[[1]][[3]][2], "[2](3:0.460,(1:0.076,2:0.076):0.384);")

  # Works on trees without recombination
  trio_trees <- generate_trio_trees(list("(2:0.865,(1:0.015,3:0.015):0.850);"),
                                    matrix(c(9, 2, 2, 5, 2, 1), 1, 6))
  expect_equal(trio_trees[[1]][[1]][1], "[9](2:0.865,(1:0.015,3:0.015):0.850);")
  expect_equal(trio_trees[[1]][[2]][1], "[2](2:0.865,(1:0.015,3:0.015):0.850);")
  expect_equal(trio_trees[[1]][[3]][1], "[2](2:0.865,(1:0.015,3:0.015):0.850);")


  # Works for non-trio loci
  trio_trees <- generate_trio_trees(trees, matrix(c(0, 0, 20, 0, 0, 2), 1, 6))
  expect_equal(length(trio_trees[[1]]$left), 0)
  expect_equal(length(trio_trees[[1]]$right), 0)
  expect_equal(length(trio_trees[[1]]$middle), 6)

  expect_error(generate_trio_trees(trees, matrix(c(0, 0, 20, 0, 0, 2), 3, 6)))
})


test_that("simulating trees for seq-gen works", {
  model <- generate_tree_model(model_theta_tau() +
                                 feat_recombination(1) +
                                 locus_averaged(2, 10) +
                                 locus_trio(locus_length = c(1,3,5),
                                            distance = c(2,4)))

  expect_true(requires_trees(model))
  stats <- simulate(model, pars = c(1, 5))
  expect_equal(length(stats$trees), 3)

  for (i in 1:2) {
    expect_equal(length(stats$trees[[i]]), 1)
    expect_true(file.exists(stats$trees[[i]]$middle))
    expect_true(file.info(stats$trees[[i]]$middle)$size > 1)
    unlink(stats$trees[[i]])
  }

  for (i in 1:3) {
    expect_true(file.exists(stats$trees[[3]][i]))
    expect_true(file.info(stats$trees[[3]][i])$size > 1)
  }
  unlink(stats$trees[[3]])
})


test_that("simulating and importing trees works", {
  model <- model_theta_tau() +
    feat_recombination(.1) +
    locus_single(10) +
    sumstat_trees("trees")

  stats <- simulate(model, pars = c(1, 5))
  expect_that(stats$trees, is_a("list"))
  expect_equal(length(stats$trees), get_locus_number(model))


  model <- model_theta_tau() + sumstat_trees("trees")
  stats <- simulate(model, pars = c(1, 5))
  expect_that(stats$trees, is_a("list"))
  expect_equal(length(stats$trees), get_locus_number(model))
  expect_true(all(sapply(stats$trees, length) == 1))
})
