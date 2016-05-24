context("BSGD")


test_that(" budget in BSGD works ", {
# convert iris to matrix
	x = as.matrix(iris[,1:4])
	y = as.vector(as.numeric(iris[,5]))
	# make sure its binary
	y = replace(y, y == 2, 0)
	y = replace(y, y == 3, 0)

	model = SharkBSGDTrain (x, y, C = 0.0001, budget = 5, gamma = 1, epochs = 1, strategy = "Merge")
	expect_equal (nrow(model$SV), 5)
})

