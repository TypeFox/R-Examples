context("Test binary")

b <- binary(byte())
bs <- binary(byte(), signed=TRUE)
bl <- binary(byte(), littleEndian=TRUE)

test_that("Lost Attributes", {
    expect_that(class(b), equals(c("binary","logical")))
    expect_that(attr(bs, "signed"), equals(TRUE))
    expect_that(attr(bl, "littleEndian"), equals(TRUE))
})

n <- c(FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE)
class(n) <- c("binary","logical")
attr(n, "signed") <- FALSE
attr(n, "littleEndian") <- FALSE
    
test_that("Return binary", {
    expect_that(length(binary(byte())), equals(byte()))
    expect_that(length(binary(byte(), signed=TRUE)), equals(byte()))
    expect_that(length(binary(1, signed=TRUE)), equals(byte()))
    expect_that(binary(byte()), equals(n))
})

context("Test as.binary")

b <- as.binary(rep(FALSE,8), logic=TRUE)
bs <- as.binary(rep(FALSE,8), signed=TRUE, logic=TRUE)
bs2 <- as.binary(rep(FALSE,1), signed=TRUE, logic=TRUE)
bl <- as.binary(rep(FALSE,8), littleEndian=TRUE, logic=TRUE)
bl2b <- binary(byte(), littleEndian=TRUE)
bl2b <- as.binary(bl2b, littleEndian=FALSE)
bs2u <- binary(byte(), signed=TRUE)
bs2u <- as.binary(bs2u, signed=FALSE)

size1BEsigned <- as.binary(2, signed=TRUE, size=1)
size2BEsigned <- as.binary(2, signed=TRUE, size=2)
size3BEsigned <- as.binary(2, signed=TRUE, size=3)
size4BEsigned <- as.binary(2, signed=TRUE, size=4)
n1BE <- as.binary(1, n=1)
n2BE <- as.binary(1, n=2)
n3BE <- as.binary(1, n=3)
n4BE <- as.binary(1, n=4)
n5BE <- as.binary(1, n=5)
n6BE <- as.binary(1, n=6)
n7BE <- as.binary(1, n=7)
n8BE <- as.binary(1, n=8)

test_that("Return size as.binary", {
    expect_that(length(size1BEsigned), equals(byte()))
    expect_that(length(size2BEsigned), equals(2*byte()))
    expect_that(length(size3BEsigned), equals(3*byte()))
    expect_that(length(size4BEsigned), equals(4*byte()))    
    expect_that(length(n1BE), equals(1))
    expect_that(length(n2BE), equals(2))
    expect_that(length(n3BE), equals(3))
    expect_that(length(n4BE), equals(4))
    expect_that(length(n5BE), equals(5))
    expect_that(length(n6BE), equals(6))
    expect_that(length(n7BE), equals(7))
    expect_that(length(n8BE), equals(8))
})

test_that("Lost Attributes", {
    expect_that(class(b), equals(c("binary","logical")))
    expect_that(attr(b, "signed"), equals(FALSE))
    expect_that(attr(b, "littleEndian"), equals(FALSE))
    expect_that(attr(bs, "signed"), equals(TRUE))
    expect_that(attr(bl, "littleEndian"), equals(TRUE))
    expect_that(attr(bl2b, "littleEndian"), equals(FALSE))
    expect_that(attr(bs2u, "signed"), equals(FALSE))    
})

test_that("Return as.binary", {
    expect_that(length(b), equals(byte()))
    expect_that(length(bs2), equals(byte()))    
})

context("Test some converting")

test_that("Return value as.binary", {
    expect_that(as.binary(0), is_a("binary"))
    expect_that(as.binary(0), is_equivalent_to(as.binary(c(0), logic=TRUE)))
    expect_that(as.binary(0, signed=TRUE, size=1), is_equivalent_to(as.binary(c(0,0,0,0,0,0,0,0), logic=TRUE)))
    expect_that(as.binary(0, littleEndian=TRUE, signed=TRUE, size=1), is_equivalent_to(as.binary(c(0,0,0,0,0,0,0,0), logic=TRUE)))
    expect_that(as.binary(1), is_equivalent_to(as.binary(c(1), signed=FALSE, littleEndian=FALSE, logic=TRUE)))
    expect_that(as.binary(1, signed=TRUE, size=1), is_equivalent_to(as.binary(c(0,0,0,0,0,0,0,1), logic=TRUE)))
    expect_that(as.binary(1, littleEndian=TRUE, signed=TRUE, size=1), is_equivalent_to(as.binary(c(1,0,0,0,0,0,0,0), signed=TRUE, littleEndian=TRUE, logic=TRUE)))
    expect_that(as.binary(-1,signed=TRUE), is_equivalent_to(as.binary(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1), logic=TRUE)))
    expect_that(as.binary(-1, littleEndian=TRUE, signed=TRUE, size=2), is_equivalent_to(as.binary(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1), logic=TRUE)))
    expect_that(as.binary(8, signed=TRUE, size=1) , is_equivalent_to(as.binary(c(0,0,0,0,1,0,0,0), signed=TRUE, logic=TRUE)))
    expect_that(as.binary(8, littleEndian=TRUE, signed=TRUE, size=1) , is_equivalent_to(as.binary(c(0,0,0,1,0,0,0,0), signed=TRUE, littleEndian=TRUE, logic=TRUE)))
})

test_that("Return value as.binary2", {
    expect_that(as.binary(0:1), is_a("list"))
    expect_that(as.binary(0:1)[[1]], is_a("binary"))
    expect_that(as.binary(0:1)[[2]], is_a("binary"))
    expect_that(length(as.binary(0:1)), equals(2))
    expect_that(as.binary(0:4)[[1]], is_equivalent_to(as.binary(c(0), logic=TRUE)))
    expect_that(as.binary(0:4)[[2]], is_equivalent_to(as.binary(c(1), logic=TRUE)))
    expect_that(as.binary(0:4)[[3]], is_equivalent_to(as.binary(c(1,0), logic=TRUE)))
    expect_that(as.binary(0:4)[[4]], is_equivalent_to(as.binary(c(1,1), logic=TRUE)))
    expect_that(as.binary(0:4)[[5]], is_equivalent_to(as.binary(c(1,0,0), logic=TRUE)))
    expect_that(as.binary(0:-4, signed=TRUE, size=1)[[1]], is_equivalent_to(as.binary(c(0,0,0,0,0,0,0,0), logic=TRUE)))
    expect_that(as.binary(0:-4, signed=TRUE, size=1)[[2]], is_equivalent_to(as.binary(c(1,1,1,1,1,1,1,1), logic=TRUE)))
    expect_that(as.binary(0:-4, signed=TRUE, size=1)[[3]], is_equivalent_to(as.binary(c(1,1,1,1,1,1,1,0), logic=TRUE)))
    expect_that(as.binary(0:-4, signed=TRUE, size=1)[[4]], is_equivalent_to(as.binary(c(1,1,1,1,1,1,0,1), logic=TRUE)))
    expect_that(as.binary(0:-4, signed=TRUE, size=1)[[5]], is_equivalent_to(as.binary(c(1,1,1,1,1,1,0,0), logic=TRUE)))
})

test_that("Return value as.numeric", {
    expect_that(as.numeric(as.binary(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1), signed=TRUE, logic=TRUE)), is_equivalent_to(-1))
    expect_that(as.numeric(as.binary(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1), signed=TRUE, littleEndian=TRUE, logic=TRUE)), is_equivalent_to(-1))
    expect_that(as.numeric(as.binary(c(0), logic=TRUE)), is_equivalent_to(0))
    expect_that(as.numeric(as.binary(c(0), littleEndian=TRUE), logic=TRUE), is_equivalent_to(0))
    expect_that(as.numeric(as.binary(c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0), signed=TRUE, littleEndian=TRUE, logic=TRUE)), is_equivalent_to(0))
    expect_that(as.numeric(as.binary(c(1), logic=TRUE)), is_equivalent_to(1))
    expect_that(as.numeric(as.binary(c(1), littleEndian=TRUE), logic=TRUE), is_equivalent_to(1))
    expect_that(as.numeric(as.binary(c(1,0,0,0,0,0,0,0), signed=TRUE, littleEndian=TRUE, logic=TRUE)), is_equivalent_to(1))
    expect_that(as.numeric(as.binary(c(0,0,0,0,0,0,0,1), signed=TRUE, logic=TRUE)), is_equivalent_to(1))
    expect_that(as.hexmode(as.numeric(as.binary(c(1,1,1,1), logic=TRUE))), is_equivalent_to(as.hexmode("f")))
})

# in work !
inputtest <- function(s)
{
    for(i in s) {
        x <- as.numeric(as.binary(i, signed=TRUE))
        if  (x != i) return(1)
    }
    return(0)
}

test_that("Range", {
    expect_that(inputtest(-64:64), equals(0))
})

context("Test is.binary")

b <- binary(1)
bd <- as.binary(1)
bl <- as.binary(TRUE, logic=TRUE)
l <- logical(1)
i <- integer(1)
d <- double(1)
n <- numeric(1)
f <- factor(1)

test_that("Return as.binary", {
    expect_that(is.binary(b), equals(TRUE))
    expect_that(is.binary(bd), equals(TRUE))
    expect_that(is.binary(bl), equals(TRUE))
    expect_that(is.binary(l), equals(FALSE))
    expect_that(is.binary(i), equals(FALSE))
    expect_that(is.binary(d), equals(FALSE))
    expect_that(is.binary(n), equals(FALSE))
    expect_that(is.binary(f), equals(FALSE))
})

context("Test print.binary")

context("Test summary.binary")

context("Test as.raw.binary")

l <- rep(TRUE, 32)
l2 <- l
l2[25] <- FALSE
l3 <- l
l3[5] <- FALSE
r <- packBits(l)
r2 <- packBits(l2)
r3 <- packBits(l3)
b <- as.binary(0xffffffff)
b2 <- as.binary(-2, signed=TRUE, size=4)
b3 <- as.binary(-2, signed=TRUE, littleEndian=TRUE, size=4)

test_that("Return as.raw.binary", {
    expect_that(length(as.raw(b)), equals(4))
    expect_that(as.raw(b), equals(r))
    expect_that(as.raw(b2), equals(r2))
    expect_that(as.raw(b3), equals(r3))
})


context("Test as.integer.binary")

zero <- as.binary(0)
zeroL <- as.binary(0, littleEndian=TRUE)
zeroS <- as.binary(0, signed=TRUE)
zeroSL <- as.binary(0, signed=TRUE, littleEndian=TRUE)
minusone <- as.binary(-1, signed=TRUE)
minusoneL <- as.binary(-1, signed=TRUE, littleEndian=TRUE)
one <- as.binary(1)
oneL <- as.binary(1, littleEndian=TRUE)
oneS <- as.binary(1, signed=TRUE)
oneSL <- as.binary(1, signed=TRUE, littleEndian=TRUE)

test_that("Type of as.integer.binary", {
    expect_that(typeof(as.integer(zero)), equals("integer"))
    expect_that(typeof(as.integer(zeroL)), equals("integer"))
    expect_that(typeof(as.integer(zeroS)), equals("integer"))
    expect_that(typeof(as.integer(zeroSL)), equals("integer"))
    expect_that(typeof(as.integer(minusone)), equals("integer"))
    expect_that(typeof(as.integer(minusoneL)), equals("integer"))
    expect_that(typeof(as.integer(one)), equals("integer"))
    expect_that(typeof(as.integer(oneL)), equals("integer"))
    expect_that(typeof(as.integer(oneS)), equals("integer"))
    expect_that(typeof(as.integer(oneSL)), equals("integer"))
})

test_that("Return as.integer.binary", {
    expect_that(as.integer(zero), equals(0))
    expect_that(as.integer(zeroL), equals(0))
    expect_that(as.integer(zeroS), equals(0))
    expect_that(as.integer(zeroSL), equals(0))
    expect_that(as.integer(minusone), equals(-1))
    expect_that(as.integer(minusoneL), equals(-1))
    expect_that(as.integer(one), equals(1))
    expect_that(as.integer(oneL), equals(1))
    expect_that(as.integer(oneS), equals(1))
    expect_that(as.integer(oneSL), equals(1))
})

context("Test as.double.binary")

zero <- as.binary(0)
zeroL <- as.binary(0, littleEndian=TRUE)
zeroS <- as.binary(0, signed=TRUE)
zeroSL <- as.binary(0, signed=TRUE, littleEndian=TRUE)
minusone <- as.binary(-1, signed=TRUE)
minusoneL <- as.binary(-1, signed=TRUE, littleEndian=TRUE)
one <- as.binary(1)
oneL <- as.binary(1, littleEndian=TRUE)
oneS <- as.binary(1, signed=TRUE)
oneSL <- as.binary(1, signed=TRUE, littleEndian=TRUE)

test_that("Type of as.double.binary", {
    expect_that(typeof(as.double(zero)), equals("double"))
    expect_that(typeof(as.double(zeroL)), equals("double"))
    expect_that(typeof(as.double(zeroS)), equals("double"))
    expect_that(typeof(as.double(zeroSL)), equals("double"))
    expect_that(typeof(as.double(minusone)), equals("double"))
    expect_that(typeof(as.double(minusoneL)), equals("double"))
    expect_that(typeof(as.double(one)), equals("double"))
    expect_that(typeof(as.double(oneL)), equals("double"))
    expect_that(typeof(as.double(oneS)), equals("double"))
    expect_that(typeof(as.double(oneSL)), equals("double"))
})

test_that("Return as.double.binary", {
    expect_that(as.double(zero), equals(0))
    expect_that(as.double(zeroL), equals(0))
    expect_that(as.double(zeroS), equals(0))
    expect_that(as.double(zeroSL), equals(0))
    expect_that(as.double(minusone), equals(-1))
    expect_that(as.double(minusoneL), equals(-1))
    expect_that(as.double(one), equals(1))
    expect_that(as.double(oneL), equals(1))
    expect_that(as.double(oneS), equals(1))
    expect_that(as.double(oneSL), equals(1))
})

context("Test Ops.binary")

context("Test '+.binary'")

mtwo <- as.binary(-2, signed=TRUE)
mone <- as.binary(-1, signed=TRUE)
zero <- as.binary(0, signed=TRUE)
one <- as.binary(1, signed=TRUE)
two <- as.binary(2, signed=TRUE)

# signed and big endian only.
test_that("Return +", {
    expect_that(zero + zero, equals(zero))
    expect_that(zero + one, equals(one))
    expect_that(one + zero, equals(one))
    expect_that(as.numeric(one + one), equals(2))
    expect_that(as.numeric(one + mone), equals(0))
    expect_that(as.numeric(mone + one), equals(0))    
    expect_that(as.numeric(zero + mone), equals(-1))
    expect_that(as.numeric(mone + zero), equals(-1))
    expect_that(as.numeric(zero + mtwo), equals(-2))
    expect_that(as.numeric(mtwo + zero), equals(-2))
    expect_that(as.numeric(one + two), equals(3))
    expect_that(as.numeric(two + one), equals(3))
    expect_that(as.numeric(two + two), equals(4))
    expect_that(as.numeric(mtwo + mtwo), equals(-4))    
})

context("Test '-.binary'")

mtwo <- as.binary(-2, signed=TRUE)
mone <- as.binary(-1, signed=TRUE)
zero <- as.binary(0, signed=TRUE)
one <- as.binary(1, signed=TRUE)
two <- as.binary(2, signed=TRUE)

# signed and big endian only.
test_that("Return -", {
    expect_that(zero - zero, equals(zero))
    expect_that(zero - one, equals(mone))
    expect_that(one - zero, equals(one))
    expect_that(as.numeric(one - one), equals(0))
    expect_that(as.numeric(one - mone), equals(2))
    expect_that(as.numeric(mone - one), equals(-2))    
    expect_that(as.numeric(zero - mone), equals(1))
    expect_that(as.numeric(mone - zero), equals(-1))
    expect_that(as.numeric(zero - mtwo), equals(2))
    expect_that(as.numeric(mtwo - zero), equals(-2))
    expect_that(as.numeric(one - two), equals(-1))
    expect_that(as.numeric(two - one), equals(1))
    expect_that(as.numeric(two - two), equals(0))
    expect_that(as.numeric(mtwo - mtwo), equals(0))    
})

context("Test '==.binary'")

s_mtwo <- as.binary(-2, signed=TRUE)
s_mone <- as.binary(-1, signed=TRUE)
s_zero <- as.binary(0, signed=TRUE)
s_one <- as.binary(1, signed=TRUE)
s_two <- as.binary(2, signed=TRUE)
zero <- as.binary(0)
one <- as.binary(1)
two <- as.binary(2)

test_that("Return -", {
    expect_that(one == two, equals(FALSE))
    expect_that(s_mtwo == s_mtwo, equals(TRUE))
    expect_that(s_mone == s_mone, equals(TRUE))
    expect_that(s_zero == s_zero, equals(TRUE))
    expect_that(s_one == s_one, equals(TRUE))
    expect_that(s_two == s_two, equals(TRUE))    
    expect_that(zero == zero, equals(TRUE))
    expect_that(one == one, equals(TRUE))
    expect_that(two == two, equals(TRUE))
    
    expect_that(switchEndianess(s_mtwo) == s_mtwo, equals(TRUE))
    expect_that(s_mone == switchEndianess(s_mone), equals(TRUE))
    expect_that(switchEndianess(s_zero )== s_zero, equals(TRUE))
    expect_that(s_one == switchEndianess(s_one), equals(TRUE))
    expect_that(switchEndianess(s_two) == s_two, equals(TRUE))
    expect_that(zero == switchEndianess(zero), equals(TRUE))
    expect_that(switchEndianess(one) == one, equals(TRUE))
    expect_that(two == switchEndianess(two), equals(TRUE))
})

context("Test '!=.binary'")

s_mtwo <- as.binary(-2, signed=TRUE)
s_mone <- as.binary(-1, signed=TRUE)
s_zero <- as.binary(0, signed=TRUE)
s_one <- as.binary(1, signed=TRUE)
s_two <- as.binary(2, signed=TRUE)
zero <- as.binary(0)
one <- as.binary(1)
two <- as.binary(2)

test_that("Return -", {
    expect_that(one != two, equals(TRUE))
    expect_that(s_mtwo != s_mtwo, equals(FALSE))
    expect_that(s_mone != s_mone, equals(FALSE))
    expect_that(s_zero != s_zero, equals(FALSE))
    expect_that(s_one != s_one, equals(FALSE))
    expect_that(s_two != s_two, equals(FALSE))
    expect_that(zero != zero, equals(FALSE))
    expect_that(one != one, equals(FALSE))
    expect_that(two != two, equals(FALSE))
    
    expect_that(switchEndianess(s_mtwo) != s_mtwo, equals(FALSE))
    expect_that(s_mone != switchEndianess(s_mone), equals(FALSE))
    expect_that(switchEndianess(s_zero )!= s_zero, equals(FALSE))
    expect_that(s_one != switchEndianess(s_one), equals(FALSE))
    expect_that(switchEndianess(s_two) != s_two, equals(FALSE))
    expect_that(zero != switchEndianess(zero), equals(FALSE))
    expect_that(switchEndianess(one) != one, equals(FALSE))
    expect_that(two != switchEndianess(two), equals(FALSE))
})

context("Test '[.binary'")

b <- binary(byte())
sb <- binary(byte(), signed=TRUE)
lb <- binary(byte(), littleEndian=TRUE)

test_that("Lost Attributes", {
    expect_that(attr(b[length(b):1], "class"), equals(c("binary","logical")))
    expect_that(attr(b[length(b):1], "signed"), equals(FALSE))
    expect_that(attr(b[length(b):1], "littleEndian"), equals(FALSE))
    expect_that(attr(sb[length(sb):1], "signed"), equals(TRUE))
    expect_that(attr(lb[length(lb):1], "littleEndian"), equals(TRUE))
})

context("Test rev")

b <- binary(byte())
sb <- binary(byte(), signed=TRUE)
lb <- binary(byte(), littleEndian=TRUE)

test_that("Lost Attributes", {
    expect_that(attr(rev(b), "class"), equals(c("binary","logical")))
    expect_that(attr(rev(b), "signed"), equals(FALSE))
    expect_that(attr(rev(b), "littleEndian"), equals(FALSE))
    expect_that(attr(rev(sb), "signed"), equals(TRUE))
    expect_that(attr(rev(lb), "littleEndian"), equals(TRUE))
})

context("Test saveAttributes")

context("Test loadAttributes")
