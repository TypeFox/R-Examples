library (lasvmR)
context ("LASVM Regressiontest")	


test_that("our wrapper gives the same output as LASVM called by command line", {
	# generate synthetic data set (see some stack overflow question)

	# generate original model via ./la_svm -o 0 -t 2 -s 1 -p 5 -g 2 -c 2 ./synthetical.sparse.train lasvm.model
	# and then predictions via ./la_test ./synthetical.sparse.test ./lasvm.model predictions
	# actually we expect here 100%, so only need to test against the test data labels.
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

	model = lasvmTrain (x = traindata, y = trainlabel, gamma = 2, cost = 2, epochs = 5, optimizer = 0, kernel = 2, selection = 1, verbose = FALSE)
	predictions = lasvmPredict (testdata, model, verbose = FALSE)
	expect_equal (sum(abs(predictions$predictions - testlabel)), 0)

	
	
	# take 0815 iris set
	set.seed(32)
	d = iris[sample(nrow(iris)),]
	x = as.matrix(d[,1:4])
	y = as.matrix(as.numeric(d[,5]))
	y[y==3] = 1
	y[y==2] = -1

	xt = as.matrix(x[101:150,])
	yt = as.matrix(as.numeric(y[101:150,]))

	model = lasvmTrain (x[1:100,], y[1:100], gamma = 0.1, cost = 1, epochs = 5, optimizer = 0, kernel = 2, selection = 1, verbose = FALSE)
	predictions = lasvmPredict (xt, model, verbose = FALSE)

	
 	# test the same for polynomial kernel
 	set.seed(32)
	degree = 3.5
	coef0 = 8.86
	cost = 0.001

	model = lasvmTrain (x[1:100,], y[1:100,], degree = degree, coef0 = coef0, cost = cost,  kernel = 1, selection = 2, verbose = FALSE)
	predictions = lasvmPredict (xt, model, verbose = FALSE)

	expect_equal (sum(abs(predictions$predictions - yt)), 0)
})


test_that("different parameter combinations for RBF kernel does not make lasvmR crash", {

	# take 0815 iris set
	d = iris[sample(nrow(iris)),]
	x = as.matrix(d[,1:4])
	y = as.matrix(as.numeric(d[,5]))
	y[y==3] = 1
	y[y==2] = -1

	for (i in 1:32) {
		s = sample(nrow(iris))
		p = round(runif(1)*(nrow(iris)-10))+5
		trIdx = s[1:p]
		testIdx = s[p+1:nrow(iris)]

		model = lasvmTrain (x[trIdx,], y[trIdx,], degree = runif(1)*10000, coef0 = runif(1)*10000, gamma = runif(1)*10000, cost = runif(1)*10000, epochs = round(runif(1)*10), optimizer = round(runif(1)), kernel = 2, selection = round(runif(1)*2), verbose = FALSE)
		predictions = lasvmPredict (x[testIdx,], model, verbose = FALSE)
	}
	
	expect_equal (1, 1)
})
