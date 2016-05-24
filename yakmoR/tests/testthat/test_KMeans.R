context("Simple Gaussian")	


# this is more like a sanity check.
test_that("KMeans trains as expected on simple 2d-gaussian dataset", {

	# generate 2 clusters
	set.seed(101)
	qx = rnorm(100, mean = -3, sd = 1) - 1
	qy = rnorm(100, mean = -3, sd = 1) - 1
	px = rnorm(100, mean = 3, sd = 1) + 1
	py = rnorm(100, mean = 3, sd = 1) + 1
	data = rbind( cbind(px, py), cbind(qx, qy) )
	E = yakmoR::orthoKMeansTrain (x = data, k = 2, rounds = 4, verbose = FALSE)
	
	# check that everything in the upper quadrant is the same cluster and same for those below
	for (i in 1:nrow(data)) {
		if (data[i,"px"] > 0) {
			expect_equal( E$cluster[[1]][i], 0)
		} else {
			expect_equal( E$cluster[[1]][i], 1)
		}
		expect_equal( E$cluster[[4]][i], 0)
	}
})


