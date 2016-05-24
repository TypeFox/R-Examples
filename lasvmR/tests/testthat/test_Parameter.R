library (lasvmR)
context("LASVM Parametertest")	


test_that("walltime works roughly as expected", {
	# generate synthetic data set (see some stack overflow question)

	# generate original model via ./la_svm -o 0 -t 2 -s 1 -p 5 -g 2 -c 2 ./synthetical.sparse.train lasvm.model
	# and then predictions via ./la_test ./synthetical.sparse.test ./lasvm.model predictions
	# NOTE: the original LASVM (v1.1) has a bug, and will not save the predictions. you just need to insert
	# a fprintf (fp, "%d\n", (int)y); at the end of the loop in the predict routine.
	
	# generate 2 clusters
	set.seed(101)
	qx = rnorm(100, mean = -3, sd = 1) - 1
	qy = rnorm(100, mean = -3, sd = 1) - 1
	px = rnorm(100, mean = 3, sd = 1) + 1
	py = rnorm(100, mean = 3, sd = 1) + 1
	traindata = rbind( cbind(px, py), cbind(qx, qy) )
	trainlabel = sign (traindata[,1])
	
	set.seed(102)
	n = 333
	qx = rnorm(n, mean = -3, sd = 1) - 1
	qy = rnorm(n, mean = -3, sd = 1) - 1
	px = rnorm(n, mean = 3, sd = 1) + 1
	py = rnorm(n, mean = 3, sd = 1) + 1
	testdata = rbind( cbind(px, py), cbind(qx, qy) )
	testlabel = sign (testdata[,1])

	mT = system.time( lasvmTrain (x = traindata, y = trainlabel, gamma = 1, cost = 1, epochs = 99999999, 
		termination = 2, 
		sample = 5, # only 5 seconds for each iteration
		kernel = 2, 
		verbose = FALSE))
	
	expect_less_than (mT[3], 10)

})



