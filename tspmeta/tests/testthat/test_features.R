context("features")

test_that("features", {
	seeds = 1:3
	sizes = c(5, 10, 100)
	for (size in sizes) {
	  for (seed in seeds) {
            set.seed(seed)
		    x = random_instance(d = 2, size = size)
		    feats = features(x)
            expect_is(feats, "list", info = paste("Features not a list for size =",
                size, "and seed =", seed))
	  }
	}
})


