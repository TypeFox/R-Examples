## ------------------------------------------------------------------------
library(recosystem)
set.seed(123) # This is a randomized algorithm
trainset = system.file("dat", "smalltrain.txt", package = "recosystem")
testset = system.file("dat", "smalltest.txt", package = "recosystem")
r = Reco()
opts = r$tune(trainset, opts = list(dim = c(10, 20, 30), lrate = c(0.1, 0.2),
                                    nthread = 1, niter = 10))
opts
r$train(trainset, opts = c(opts$min, nthread = 1, niter = 20))
outfile = tempfile()
r$predict(testset, outfile)

## Compare the first few true values of testing data
## with predicted ones
# True values
print(read.table(testset, header = FALSE, sep = " ", nrows = 10)$V3)
# Predicted values
print(scan(outfile, n = 10))

