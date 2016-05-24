library("R.filesets")

# Tests adapted from the IRanges package and http://www.bioconductor.org/help/course-materials/2012/SeattleFeb2012/GenomicRanges_slides.pdf

# Example files
path <- system.file("exData", "dataSetA,original", package="R.filesets")
print(path)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setting up a file set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ds <- GenericDataFileSet$byPath(path)
print(ds)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Vector operations on data set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Single-bracket subsetting
dsHead <- ds[1:2]
dsTail <- ds[-(1:2)]

# Combining: c()
ds2 <- c(dsHead, dsTail)
stopifnot(equals(ds2, ds))



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# List operations on data set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Double-bracket subsetting
dfA <- ds[["1.2(a)"]]
dfB <- ds[["fileD,3cols"]]

# equals()
stopifnot(equals(dfA, dfA))
stopifnot(equals(dfB, dfB))
stopifnot(!equals(dfA, dfB))

# Comparing: ==, !=
# FIXME

# Ordering: <=, >=, <, >, order(), sort(), rank()
# FIXME


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Vector and list operations on data set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dsA <- extract(ds, c(seq_along(ds), 1, 3))
stopifnot(length(dsA) == length(ds) + 2L)

dsB <- c(ds, ds[[1]], ds[[3]])
stopifnot(equals(dsB, dsA))

dsC <- c(ds, list(ds[[1]], ds[[3]]))
stopifnot(equals(dsC, dsA))

dsD <- c(ds, ds[c(1,3)])
stopifnot(equals(dsD, dsA))

# Comparing: duplicated(), where duplication test is based on equals()
dups <- duplicated(dsA)
print(dups)
hasDups <- anyDuplicated(dsA)
print(hasDups)
stopifnot(identical(any(dups), hasDups))

dsT <- dsA[!dups]
stopifnot(!anyDuplicated(dsT))

# Comparing: unique()
dsU <- unique(dsA)
stopifnot(equals(dsU, dsT))
