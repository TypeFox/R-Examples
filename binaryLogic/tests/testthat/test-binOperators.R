context("Test negate")

test_that("Lost Attributes", {
    expect_that(attr(negate(binary(byte(), signed=TRUE)), "class"), equals(c("binary","logical")))
    expect_that(attr(negate(as.binary(2, signed=TRUE)), "littleEndian"), equals(FALSE))
    expect_that(attr(negate(as.binary(2, signed=TRUE)), "signed"), equals(TRUE))
    expect_that(attr(negate(as.binary(2, signed=TRUE, littleEndian=TRUE)), "littleEndian"), equals(TRUE))
})

test_that("Return negate", {
    expect_that(as.numeric(negate(binary(byte(), signed=TRUE))), equals(0))
    expect_that(as.numeric(negate(as.binary(0, signed=TRUE))), equals(0))
    expect_that(as.numeric(negate(as.binary(0, signed=FALSE))), equals(0))
    expect_that(as.numeric(negate(as.binary(0, signed=TRUE, littleEndian=TRUE))), equals(0))
    expect_that(as.numeric(negate(as.binary(0, signed=FALSE, littleEndian=TRUE))), equals(0))
    expect_that(as.numeric(negate(as.binary(-1, signed=TRUE))), equals(1))
    expect_that(as.numeric(negate(as.binary(-1, signed=TRUE, littleEndian=TRUE))), equals(1))
    expect_that(as.numeric(negate(as.binary(1, signed=TRUE))), equals(-1))
    expect_that(as.numeric(negate(as.binary(1, signed=FALSE))), equals(-1))
    expect_that(as.numeric(negate(as.binary(1, signed=TRUE, littleEndian=TRUE))), equals(-1))
    expect_that(as.numeric(negate(as.binary(1, signed=FALSE, littleEndian=TRUE))), equals(-1))
})

context("Test shiftLeft")

test_that("Lost Attributes", {
    expect_that(class(shiftLeft(binary(byte()),1)), equals(c("binary","logical")))
    expect_that(attr(shiftLeft(binary(byte()),1), "signed"), equals(FALSE))
    expect_that(attr(shiftLeft(binary(byte()),1), "littleEndian"), equals(FALSE))
    expect_that(attr(shiftLeft(binary(byte(), signed=TRUE),1), "signed"), equals(TRUE))
    expect_that(attr(shiftLeft(binary(byte(), littleEndian=TRUE),1), "littleEndian"), equals(TRUE))
    expect_that(class(shiftLeft(logical(byte()),1)), equals("logical"))
})

one <- as.binary(1, signed=TRUE, size=1)
l <- rep(FALSE, 4)
l2 <- c(TRUE,TRUE,FALSE,TRUE)
l3 <- c(TRUE,FALSE,TRUE,FALSE)
l4 <- c(FALSE,TRUE,FALSE,FALSE)
l5 <- c(TRUE,FALSE,FALSE,FALSE)

test_that("Return shiftLeft", {
    expect_that(length(shiftLeft(binary(byte()), 1)), equals(byte()))
    expect_that(length(shiftLeft(logical(byte()), 1)), equals(byte()))
    expect_that(as.numeric(shiftLeft(one, 1)), equals(2))
    expect_that(as.numeric(shiftLeft(one, 2)), equals(4))
    expect_that(as.numeric(shiftLeft(one, 3)), equals(8))
    expect_that(shiftLeft(l, 1), equals(l))
    expect_that(shiftLeft(l2, 1), equals(l3))
    expect_that(shiftLeft(l2, 2), equals(l4))
    expect_that(shiftLeft(l2, 3), equals(l5))
    expect_that(class(shiftLeft(logical(byte()),byte() + 1)), equals("logical"))
})

context("Test shiftRight")

test_that("Lost Attributes", {
    expect_that(class(shiftRight(binary(byte()),1)), equals(c("binary","logical")))
    expect_that(attr(shiftRight(binary(byte()),1), "signed"), equals(FALSE))
    expect_that(attr(shiftRight(binary(byte()),1), "littleEndian"), equals(FALSE))
    expect_that(attr(shiftRight(binary(byte(), signed=TRUE),1), "signed"), equals(TRUE))
    expect_that(attr(shiftRight(binary(byte(), littleEndian=TRUE),1), "littleEndian"), equals(TRUE))
    expect_that(class(shiftRight(logical(byte()),1)), equals("logical"))
})

eight <- as.binary(8, signed=TRUE, size=1)
l <- rep(FALSE, 4)
l2 <- c(TRUE,TRUE,FALSE,TRUE)
l3 <- c(FALSE,TRUE,TRUE,FALSE)
l4 <- c(FALSE,FALSE,TRUE,TRUE)
l5 <- c(FALSE,FALSE,FALSE,TRUE)

test_that("Return shiftRight", {
    expect_that(length(shiftRight(binary(byte()), 1)), equals(byte()))
    expect_that(as.numeric(shiftRight(eight, 1)), equals(4))
    expect_that(as.numeric(shiftRight(eight, 2)), equals(2))
    expect_that(as.numeric(shiftRight(eight, 3)), equals(1))
    expect_that(as.numeric(shiftRight(eight, 4)), equals(0))
    expect_that(shiftRight(l, 1), equals(l))
    expect_that(shiftRight(l2, 1), equals(l3))
    expect_that(shiftRight(l2, 2), equals(l4))
    expect_that(shiftRight(l2, 3), equals(l5))
    expect_that(class(shiftRight(logical(byte()),byte() + 1)), equals("logical"))
})

context("Test rotate")

test_that("Lost Attributes", {
    expect_that(class(rotate(binary(byte()),1)), equals(c("binary","logical")))
    expect_that(attr(rotate(binary(byte()),1), "signed"), equals(FALSE))
    expect_that(attr(rotate(binary(byte()),1), "littleEndian"), equals(FALSE))
    expect_that(attr(rotate(binary(byte(), signed=TRUE),1), "signed"), equals(TRUE))
    expect_that(attr(rotate(binary(byte(), littleEndian=TRUE),1), "littleEndian"), equals(TRUE))
    expect_that(class(rotate(logical(byte()),1)), equals("logical"))
})

one <- as.binary(1, signed=TRUE, size=1)
l <- rep(FALSE, 4)
l2 <- c(TRUE,TRUE,FALSE,TRUE)
l3 <- c(TRUE,FALSE,TRUE,TRUE)
l4 <- c(FALSE,TRUE,TRUE,TRUE)
l5 <- c(TRUE,TRUE,TRUE,FALSE)

test_that("Return rotate", {
    expect_that(length(rotate(binary(byte()), 1)), equals(byte()))
    expect_that(as.numeric(rotate(one, 1)), equals(2))
    expect_that(as.numeric(rotate(one, 2)), equals(4))
    expect_that(as.numeric(rotate(one, 3)), equals(8))
    expect_that(rotate(l, 1), equals(l))
    expect_that(rotate(l2, 1), equals(l3))
    expect_that(rotate(l2, 2), equals(l4))
    expect_that(rotate(l2, 3), equals(l5))      
})

context("Test fillUpToByte")

test_that("Lost Attributes", {
    expect_that(class(fillUpToByte(binary(byte()), 2)), equals(c("binary","logical")))
    expect_that(attr(fillUpToByte(binary(byte()), 2), "signed"), equals(FALSE))
    expect_that(attr(fillUpToByte(binary(byte()), 2), "littleEndian"), equals(FALSE))
    expect_that(attr(fillUpToByte(binary(byte(), signed=TRUE), 2), "signed"), equals(TRUE))
    expect_that(attr(fillUpToByte(binary(byte(), littleEndian=TRUE), 2), "littleEndian"), equals(TRUE))
})

input1 <- as.binary(c(TRUE,TRUE,FALSE,TRUE), logic=TRUE)
l2 <- as.binary(c(FALSE,FALSE,FALSE,FALSE,TRUE,TRUE,FALSE,TRUE), logic=TRUE)
l3 <- as.binary(c(TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,FALSE,TRUE), logic=TRUE)
input2 <- as.binary(c(TRUE,TRUE,FALSE,TRUE), littleEndian=TRUE, logic=TRUE)
l4 <- as.binary(c(TRUE,TRUE,FALSE,TRUE,FALSE,FALSE,FALSE,FALSE), littleEndian=TRUE, logic=TRUE)
l5 <- as.binary(c(TRUE,TRUE,FALSE,TRUE,TRUE,TRUE,TRUE,TRUE), littleEndian=TRUE, logic=TRUE)

test_that("Return fillUpToByte", {
    expect_that(fillUpToByte(input1, value=FALSE, size=1), equals(l2))
    expect_that(fillUpToByte(input1, value=TRUE, size=1), equals(l3))
    expect_that(fillUpToByte(input2, value=FALSE, size=1), equals(l4))
    expect_that(fillUpToByte(input2, value=TRUE, size=1), equals(l5))    
})

context("Test switchEndianess")

test_that("Lost Attributes", {
    expect_that(class(switchEndianess(binary(byte()))), equals(c("binary","logical")))
    expect_that(attr(switchEndianess(binary(byte())), "signed"), equals(FALSE))
    expect_that(attr(switchEndianess(binary(byte())), "littleEndian"), equals(TRUE))
    expect_that(attr(switchEndianess(binary(byte(), signed=TRUE)), "signed"), equals(TRUE))
    expect_that(attr(switchEndianess(binary(byte(), littleEndian=TRUE)), "littleEndian"), equals(FALSE))
})

s1 <- as.binary(c(FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,TRUE,TRUE), logic=TRUE)
s2 <- as.binary(c(TRUE,TRUE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE), littleEndian=TRUE, logic=TRUE)
s3 <- as.binary(c(FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,TRUE,TRUE), signed=TRUE, logic=TRUE)
s4 <- as.binary(c(TRUE,TRUE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE), signed=TRUE, littleEndian=TRUE, logic=TRUE)
s5 <- as.binary(c(FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,TRUE,TRUE), logic=TRUE)
s6 <- as.binary(c(FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,TRUE,TRUE), littleEndian=TRUE, logic=TRUE)

test_that("Return switchEndianess", {
    expect_that(switchEndianess(s1), equals(s2))
    expect_that(switchEndianess(s2), equals(s1))
    expect_that(switchEndianess(s3), equals(s4))
    expect_that(switchEndianess(s4), equals(s3))
    expect_that(switchEndianess(s5, stickyBits=TRUE), equals(s6))
    expect_that(switchEndianess(s6, stickyBits=TRUE), equals(s5))
})

context("Test binarySeq")
