require("sfsmisc")

options(warn=2)
AsciiToInt(LETTERS) # gave '.. embedded nul ..' warning

## just for fun -- typically shows "iso-latin1 charset
cat(chars8bit(1:255),"\"\n")

## Checking the new 'ndigits' default argument for digitsBase():
ee <- 0:30
for(base in 2:64)
    stopifnot((be <- base^ee) > 0, any(ok <- be < 2^52),
	      ee == floor(1e-9+ log(be, base)),
	      be[ok] == as.integer(digitsBase(be[ok], base=base)))
## failed, e.g. for 3^5, in sfsmisc <= 1.0-22

## Tests for is.whole (taken from the examples)
set.seed(307)
a <- array(runif(24), dim = c(2, 3, 4))
a[4:8] <- 4:8
m <- matrix(runif(12), 3, 4)
m[2:4] <- 2:4
v <- complex(real = seq(0.5, 1.5, by = 0.1),
             imaginary = seq(2.5, 3.5, by = 0.1))

## Find whole entries
stopifnot(identical(is.whole(a), a == round(a)),
          identical(is.whole(m), m == round(m)),
          which(is.whole(v)) == 6)

## Numbers of class integer are always whole
stopifnot(is.whole(dim(a)), is.whole(length(v)), is.whole(-1L))
