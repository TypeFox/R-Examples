library("R.filesets")

message("*** readDataFrame()")

path <- system.file("exData", "dataSetA,original", package="R.filesets")
pathnames <- list.files(path=path, pattern="[.]txt$", full.names=TRUE)
pathname <- pathnames[1]

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Basic reading
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
data <- readDataFrame(pathname)
print(data)

data <- readDataFrame(basename(pathname), path=dirname(pathname))
print(data)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Reading gzip'ed file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
pathT <- tempdir()
pathnameZ <- file.path(pathT, sprintf("%s.gz", basename(pathname)))
R.utils::gzip(pathname, pathnameZ, remove=FALSE)
dataZ <- readDataFrame(pathnameZ)
print(dataZ)

## Validate
stopifnot(identical(dataZ, data))

## Cleanup
file.remove(pathnameZ)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Reading multiple files and stack them
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
pathnames <- rep(pathname, times=3L)
data <- readDataFrame(pathnames)
print(data)
