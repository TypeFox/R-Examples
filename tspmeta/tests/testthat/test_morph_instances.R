context("morph_instances")

test_that("morph_instances", {
  seeds = 1:10
  alphas = c(0, 0.3, 1)
  size = 10
  for (alpha in alphas) {
    for (seed in seeds) {
      set.seed(seed)
      x = random_instance(size = size)
      y = random_instance(size = size)
      z = morph_instances(x, y, alpha = alpha)
      expect_is(z, "tsp_instance", info = paste("morphed element is no tsp_instance object for alpha =",
        alpha, "and seed =", seed))
    }
  }
})
