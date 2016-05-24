context("powerbydesign functions")

test_that("design.anova checks", {
  expect_that(design.anova(), throws_error())
})

test_that("boot.power.anova checks", {
  # edit if design.anova function changes!
  # manually define design
  design <- list(
    between = list(age = c("young","old")),
    within = list(condition = c("cond1","cond2","cond3")),
    num_between_conds = 2,
    factor_names = c("condition","age"),
    cov_matrix = matrix(c(4.0,1.4,1.4,0.0,0.0,0.0,1.4,1.0,0.7,0.0,0.0,0.0,1.4,0.7,1.0,0.0,0.0,0.0,0.0,0.0,0.0,4.0,2.8,1.4,0.0,0.0,0.0,2.8,4.0,1.4,0.0,0.0,0.0,1.4,1.4,1.0),nrow=6,
                        dimnames = list(c("cond1.young","cond2.young","cond3.young","cond1.old","cond2.old","cond3.old"),c("cond1.young","cond2.young","cond3.young","cond1.old","cond2.old","cond3.old"))),
    means = c(1,2,3,1,2,8)
  )
  class(design) <- "design.anova"

  expect_that(boot.power.anova(), throws_error())
  expect_that(boot.power.anova(design), throws_error())
  expect_that(boot.power.anova(design,n_from=2), throws_error())
  expect_that(boot.power.anova(design,n_from=2,n_to=4), throws_error())
  expect_that(boot.power.anova(design,n_from=2,n_to=4,num_iterations_bootstrap=10), gives_warning())
  expect_that(boot.power.anova(design,n_from=2,n_to=4,num_iterations_bootstrap=10), is_a("power_by_samplesize"))
})

