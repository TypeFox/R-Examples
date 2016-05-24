context("Test binaryPrefix")

test_that("Return value", {
    expect_that(binaryPrefix(0.5, "KiB"), is_equivalent_to(512))
    expect_that(binaryPrefix(1, "KiB"), is_equivalent_to(2^10))
    expect_that(binaryPrefix(1, "MiB"), is_equivalent_to(2^20))
    expect_that(binaryPrefix(1, "GiB"), is_equivalent_to(2^30))
    expect_that(binaryPrefix(1, "TiB"), is_equivalent_to(2^40))
    expect_that(binaryPrefix(1, "PiB"), is_equivalent_to(2^50))
    expect_that(binaryPrefix(1, "EiB"), is_equivalent_to(2^60))
    expect_that(binaryPrefix(1, "ZiB"), is_equivalent_to(2^70))
    expect_that(binaryPrefix(1, "YiB"), is_equivalent_to(2^80))
})

context("Test bytesNeeded")

context("Test byte")

test_that("Return value", {
    expect_that(byte(), is_identical_to(8))
})

context("Text gray code")

test_that("Return value bin2gray | logial vector | always bigEndian", {
    expect_that(bin2gray(as.binary(0)), equals(c(FALSE)))
    expect_that(bin2gray(as.binary(1)), equals(c(TRUE)))
    expect_that(bin2gray(as.binary(0, littleEndian=TRUE)), equals(c(FALSE)))
    expect_that(bin2gray(as.binary(1, littleEndian=TRUE)), equals(c(TRUE)))  
})

test_that("Return value bin2gray | logial vector | always bigEndian", {
    expect_that(bin2gray(as.logical(c(0,0,0))), equals(c(FALSE,FALSE,FALSE)))
    expect_that(bin2gray(as.logical(c(0,0,1))), equals(c(FALSE,FALSE,TRUE)))
    expect_that(bin2gray(as.logical(c(0,1,0))), equals(c(FALSE,TRUE,TRUE)))
    expect_that(bin2gray(as.logical(c(0,1,1))), equals(c(FALSE,TRUE,FALSE)))
    expect_that(bin2gray(as.logical(c(1,0,0))), equals(c(TRUE,TRUE,FALSE)))
    expect_that(bin2gray(as.logical(c(1,0,1))), equals(c(TRUE,TRUE,TRUE)))
    expect_that(bin2gray(as.logical(c(1,1,0))), equals(c(TRUE,FALSE,TRUE)))
    expect_that(bin2gray(as.logical(c(1,1,1))), equals(c(TRUE,FALSE,FALSE)))    
})

test_that("Return value bin2gray | binary number | bigEndian", {
    expect_that(bin2gray(fillUpToBit(as.binary(0), 3)), equals(c(FALSE,FALSE,FALSE)))
    expect_that(bin2gray(fillUpToBit(as.binary(1), 3)), equals(c(FALSE,FALSE,TRUE)))
    expect_that(bin2gray(fillUpToBit(as.binary(2), 3)), equals(c(FALSE,TRUE,TRUE)))
    expect_that(bin2gray(fillUpToBit(as.binary(3), 3)), equals(c(FALSE,TRUE,FALSE)))
    expect_that(bin2gray(fillUpToBit(as.binary(4), 3)), equals(c(TRUE,TRUE,FALSE)))
    expect_that(bin2gray(fillUpToBit(as.binary(5), 3)), equals(c(TRUE,TRUE,TRUE)))
    expect_that(bin2gray(fillUpToBit(as.binary(6), 3)), equals(c(TRUE,FALSE,TRUE)))
    expect_that(bin2gray(fillUpToBit(as.binary(7), 3)), equals(c(TRUE,FALSE,FALSE)))    
})

test_that("Return value bin2gray | littleEndian | littleEndian", {
    expect_that(bin2gray(fillUpToBit(as.binary(0, littleEndian=TRUE), 3)), equals(c(FALSE,FALSE,FALSE)))
    expect_that(bin2gray(fillUpToBit(as.binary(1, littleEndian=TRUE), 3)), equals(c(TRUE,FALSE,FALSE)))
    expect_that(bin2gray(fillUpToBit(as.binary(2, littleEndian=TRUE), 3)), equals(c(TRUE,TRUE,FALSE)))
    expect_that(bin2gray(fillUpToBit(as.binary(3, littleEndian=TRUE), 3)), equals(c(FALSE,TRUE,FALSE)))
    expect_that(bin2gray(fillUpToBit(as.binary(4, littleEndian=TRUE), 3)), equals(c(FALSE,TRUE,TRUE)))
    expect_that(bin2gray(fillUpToBit(as.binary(5, littleEndian=TRUE), 3)), equals(c(TRUE,TRUE,TRUE)))
    expect_that(bin2gray(fillUpToBit(as.binary(6, littleEndian=TRUE), 3)), equals(c(TRUE,FALSE,TRUE)))
    expect_that(bin2gray(fillUpToBit(as.binary(7, littleEndian=TRUE), 3)), equals(c(FALSE,FALSE,TRUE)))    
})

b <- as.binary(0:7, n=3)
g <- lapply(b, bin2gray)
r <- lapply(g, gray2bin)

test_that("Return value gray2bin | binary number | bigEndian", {
    expect_that(r[[1]], equals(b[[1]]))
    expect_that(r[[2]], equals(b[[2]]))
    expect_that(r[[3]], equals(b[[3]]))
    expect_that(r[[4]], equals(b[[4]]))
    expect_that(r[[5]], equals(b[[5]]))
    expect_that(r[[6]], equals(b[[6]]))
    expect_that(r[[7]], equals(b[[7]]))
    expect_that(r[[8]], equals(b[[8]]))
})

b <- as.binary(0:7, littleEndian=TRUE, n=3)
g <- lapply(b, bin2gray)
r <- lapply(g, gray2bin, littleEndian=TRUE)

test_that("Return value gray2bin | binary number | bigEndian", {
    expect_that(r[[1]], equals(b[[1]]))
    expect_that(r[[2]], equals(b[[2]]))
    expect_that(r[[3]], equals(b[[3]]))
    expect_that(r[[4]], equals(b[[4]]))
    expect_that(r[[5]], equals(b[[5]]))
    expect_that(r[[6]], equals(b[[6]]))
    expect_that(r[[7]], equals(b[[7]]))
    expect_that(r[[8]], equals(b[[8]]))
})

z <- as.binary(0)
lez <- as.binary(0, littleEndian=TRUE)
o <- as.binary(1)
leo <- as.binary(1, littleEndian=TRUE)

gz <- bin2gray(z)
glez <- bin2gray(lez)
go <- bin2gray(o)
gleo <- bin2gray(leo)

rz <- gray2bin(gz)
rlez <- gray2bin(glez, littleEndian=TRUE)
ro <- gray2bin(go)
rleo <- gray2bin(gleo, littleEndian=TRUE)

test_that("Return value gray2bin | binary number | bigEndian", {
    expect_that(length(rz), equals(1))
    expect_that(rz, equals(z))
    expect_that(length(rlez), equals(1))
    expect_that(rlez, equals(lez))
    expect_that(length(ro), equals(1))
    expect_that(ro, equals(o))
    expect_that(length(rleo), equals(1))
    expect_that(rleo, equals(leo))
})