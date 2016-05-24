library("dChipIO")

path <- system.file("exData", package="dChipIO")
filename <- "Test3-1-121502.dcp"
pathname <- file.path(path, filename)

hdr <- readDcpHeader(pathname)
print(hdr)

data <- readDcp(pathname)
str(data)

# Read a subset of the units
units <- c(10:11, 15:20, 150:105, 2,2,2)
dataT <- readDcp(pathname, units=units)
str(dataT)

# Assert correctness
for (ff in c("calls", "thetas", "thetaStds", "excludes")) {
  stopifnot(length(dataT[[ff]]) == length(units))
  stopifnot(identical(dataT[[ff]], data[[ff]][units]))
}
