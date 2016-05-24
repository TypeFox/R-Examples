library(retistruct)

## Test reading file with one row in datapoints.csv file
dataset <- file.path(system.file(package = "retistruct"), "extdata", "ijroi1")
r <- retistruct.read.dataset(dataset)
print(r$Ds)

## Test reading file with two rows in datapoints.csv file
dataset <- file.path(system.file(package = "retistruct"), "extdata", "ijroi2")
r <- retistruct.read.dataset(dataset)
print(r$Ds)
