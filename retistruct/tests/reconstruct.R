library(retistruct)
dataset <- file.path(system.file(package = "retistruct"), "extdata", "GMB530/R-CONTRA")
r <- retistruct.read.dataset(dataset)
## Load the human annotation of tears
r <- retistruct.read.markup(r)
## Reconstruct
r <- retistruct.reconstruct(r)

## Save as matlab
r$dataset <- tempdir()
retistruct.export.matlab(r)
print(paste("Output saved to ", r$dataset))

