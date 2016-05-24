library(jpeg)

## grayscale
s0 <- matrix(0:9999/9999, 100)
j0 <- writeJPEG(s0, raw())
i0 <- readJPEG(j0)

# allow 2% tolerance when comparing uncompressed and compressed images
# since JPEG is lossy (the default quality is 0.7 which should be good enough)
tolerance <- 0.02

stopifnot(identical(dim(i0), dim(s0)))
# JPEG is lossy so there will be differences but they should not be too big
stopifnot(max(abs(s0 - i0)) < tolerance)

n0 <- readJPEG(j0, native=TRUE)
stopifnot(identical(dim(i0), dim(s0)))
stopifnot(inherits(n0, "nativeRaster") && identical(attr(n0, "channels"), 1L))

# check the native result for sanity - it should be XXXA
# the 8 MSB must be 1 since the alpha is 1.0 (-16777216L .. 0L)
stopifnot(all(n0 < 0L & n0 >= -16777216L))
# remove the MSB
y <- n0 + 16777216L
x <- as.integer(s0 * 255 + 0.5)
stopifnot(max(abs(x - t(y %% 256L))) < tolerance * 255)
stopifnot(all(as.integer(y / 256L) %% 256L == y %% 256L))
stopifnot(all(as.integer(y / 65536L) %% 256L == y %% 256L))

# check file vs in-memory
writeJPEG(s0, "image0.jpeg")
s <- file.info("image0.jpeg")$size
stopifnot(all(s == length(j0)))
f <- file("image0.jpeg", "rb")
j0f <- readBin(f, raw(), s)
close(f)
stopifnot(identical(c(j0f), c(j0)))
i0f <- readJPEG("image0.jpeg")
stopifnot(identical(i0f, i0))
n0f <- readJPEG("image0.jpeg", native=TRUE)
stopifnot(identical(n0f, n0))


## GA + alpha mixing
a1 <- array(c(s0, rev(s0)), c(100L, 100L, 2L))
j1 <- writeJPEG(a1, raw(), bg="black")
i1 <- readJPEG(j1)
# since JPEG flattens alpha it will have the dimensions of s0 instead of a1
stopifnot(identical(dim(i1), dim(s0)))
s1 <- s0 * rev(s0) ## this should be the result of alpha blending with black
stopifnot(max(abs(s1 - i1)) < tolerance)
i1.1 <- readJPEG(writeJPEG(a1, raw(), bg="white"))
s1.1 <- s0 * rev(s0) + (1 - rev(s0))
stopifnot(max(abs(s1.1 - i1.1)) < tolerance)


## RGB
a2 <- array(c(s0, t(s0), rev(s0)), c(100L, 100L, 3L))
j2 <- writeJPEG(a2, raw())
i2 <- readJPEG(j2)
stopifnot(identical(dim(a2), dim(i2)))
# more tolerance since we have 3x more data to compress
stopifnot(max(abs(a2 - i2)) < tolerance * 3)

# since RGB is most frequently used, check file vs raw as well
writeJPEG(a2, "image2.jpeg")
s <- file.info("image2.jpeg")$size
stopifnot(all(s == length(j2)))
f <- file("image2.jpeg", "rb")
j2f <- readBin(f, raw(), s)
close(f)
stopifnot(identical(c(j2f), c(j2)))
i2f <- readJPEG("image2.jpeg")
stopifnot(identical(i2f, i2))
n2f <- readJPEG("image2.jpeg", native=TRUE)
n2 <- readJPEG(j2, native=TRUE)
stopifnot(identical(n2f, n2))


## RGB + alpha mixing
a3 <- array(c(s0, t(s0), rev(s0), t(rev(s0))), c(100L, 100L, 4L))
j3 <- writeJPEG(a3, raw(), bg="black")
i3 <- readJPEG(j3)
# we use a2 to compare to we just added alpha
stopifnot(max(abs(i3 - a2 * rev(s0))) < tolerance * 3)
j3.1 <- writeJPEG(a3, raw(), bg="white")
i3.1 <- readJPEG(j3.1)
stopifnot(max(abs(i3.1 - a2 * rev(s0) - (1 - rev(s0)))) < tolerance * 3)

## external file checks
## those are already used in examples so it's not really necessary ..
fn <- system.file("img", "Rlogo.jpg", package="jpeg")
i4 <- readJPEG(fn)
s <- file.info(fn)$size
f <- file(fn, "rb")
j4 <- readBin(f, raw(), s)
close(f)
i4.1 <- readJPEG(fn)
stopifnot(identical(i4, i4.1))


## large RGB check
s5 <- matrix(0:999999/999999, 1000)
a5 <- array(c(s5, t(s5), rev(s5)), c(1000L, 1000L, 3L))
# produce larger files
j5 <- writeJPEG(a5, raw(), quality=0.9)
writeJPEG(a5, "image5.jpeg", quality=0.9)
s <- file.info("image5.jpeg")$size
stopifnot(all(s == length(j5)))
f <- file("image5.jpeg", "rb")
j5f <- readBin(f, raw(), s)
close(f)
stopifnot(identical(c(j5f), c(j5)))
i5 <- readJPEG(j5, native=TRUE)
i5f <- readJPEG("image5.jpeg", native=TRUE)
stopifnot(identical(i5, i5f))


## Wohoo! all tests passed!
