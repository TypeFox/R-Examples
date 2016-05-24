library("dChipIO")

path <- system.file("exData", package="dChipIO")
chipType <- "Test3"
filename <- sprintf("%s.CDF.bin", chipType)
pathname <- file.path(path, filename)

hdr <- readCdfBinHeader(pathname)
print(hdr)

data <- readCdfBin(pathname)
str(data)

# Read a subset of the units
units <- c(10:11, 15:20, 150:105, 2,2,2)
dataT <- readCdfBin(pathname, units=units)
str(dataT)

# Assert correctness
for (ff in c("unitNames", "numProbes", "CellPos")) {
  stopifnot(length(dataT[[ff]]) == length(units))
  stopifnot(identical(dataT[[ff]], data[[ff]][units]))
}
