#!/usr/bin/Rscript

library(RcppCNPy)

set.seed(42)
M1 <- matrix(rnorm(1e6), 1e3, 1e3)

## double precision floating point, uncompressed
tempmatfile <- tempfile(pattern="npymat", fileext=".npy")
npySave(tempmatfile, M1)
M2 <- npyLoad(tempmatfile)
identical(M1, M2)

## double precision floating point, compressed
tempmatfile <- tempfile(pattern="npymat", fileext=".npy.gz")
npySave(tempmatfile, M1)
M3 <- npyLoad(tempmatfile)
identical(M1, M3)

## double precision floating point vector
dim(M1) <- NULL
tempmatfile <- tempfile(pattern="npyvec", fileext=".npy")
npySave(tempmatfile, M1)
M2 <- npyLoad(tempmatfile)
identical(M1, M2)

## double precision floating point vector, compressed
tempmatfile <- tempfile(pattern="npyvec", fileext=".npy.gz")
npySave(tempmatfile, M1)
M3 <- npyLoad(tempmatfile)
identical(M1, M3)

## integer, uncompressed
M4 <- matrix(as.integer(round(M1)), 1e3, 1e3)
tempmatfile <- tempfile(pattern="intnpymat", fileext=".npy")
npySave(tempmatfile, M4)
M5 <- npyLoad(tempmatfile, "integer")
identical(M4, M5)

## integer, compressed
tempmatfile <- tempfile(pattern="intnpymat", fileext=".npy.gz")
npySave(tempmatfile, M4)
M6 <- npyLoad(tempmatfile, "integer")
identical(M4, M6)

## integer vector, uncompressed
dim(M4) <- NULL
tempmatfile <- tempfile(pattern="intnpyvec", fileext=".npy")
npySave(tempmatfile, M4)
M5 <- npyLoad(tempmatfile, "integer")
identical(M4, M5)

## integer vector, compressed
tempmatfile <- tempfile(pattern="intnpyvec", fileext=".npy.gz")
npySave(tempmatfile, M4)
M6 <- npyLoad(tempmatfile, "integer")
identical(M4, M6)
