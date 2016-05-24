library("R.filesets")
library("R.utils")  # writeDataFrame()

message("*** A huge TabularTextFile")

# Generate a large data frame
set.seed(42)
n <- 1e5
data <- data.frame(
  index = seq_len(n),
  x = runif(n),
  y = rnorm(n),
  symbol = sample(letters, size=n, replace=TRUE)
)

# Write to tab-delimited file
pathname <- tempfile(fileext=".tsv")
writeDataFrame(data, file=pathname)

# Setup tabular file
db <- TabularTextFile(pathname)
print(db)

# Read subset of the rows
n <- nbrOfRows(db)
data2 <- readDataFrame(db, rows=seq(from=1, to=n, length.out=0.10*n))
str(data2)

# Cleanup
file.remove(pathname)
