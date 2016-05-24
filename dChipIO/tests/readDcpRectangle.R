library("dChipIO")

path <- system.file("exData", package="dChipIO")
filename <- "Test3-1-121502.dcp"
pathname <- file.path(path, filename)

data <- readDcpRectangle(pathname)

layout(matrix(1:4, nrow=2, byrow=TRUE))
image(data$rawIntensities, main="Raw probe signals")
image(data$normalizedIntensities, main="Normalized probe signals")


