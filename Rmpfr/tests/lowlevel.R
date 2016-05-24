#### Low level stuff - debugging etc
#### =========         =========

require("Rmpfr")
options(warn = 2)# warning -> error

identical3 <- function(x,y,z)	  identical(x,y) && identical (y,z)
identical4 <- function(a,b,c,d)   identical(a,b) && identical3(b,c,d)

###----- _1_ mpfr1 , import, xport etc -----------------------------------------
i8 <- mpfr(-2:5, 32)
x4 <- mpfr(c(NA, NaN, -Inf, Inf), 32); x4 # NA -> NaN as well
stopifnot(identical3(is.na(x4), is.nan(x4), c(T,T,F,F)))

o1 <- as(x4[1], "mpfr1")
stopifnot(is(o1, "mpfr1")) # failed previously
validObject(o1)            # ditto (failed on 64-bit only)

stopifnot(
    getPrec("0xabc", base=16, doNumeric=FALSE) == 3*4,
    getPrec(  "abc", base=16, doNumeric=FALSE) == 3*4,
    getPrec("0b1001", base=2, doNumeric=FALSE) == 4,
    getPrec(  "1001", base=2, doNumeric=FALSE) == 4,
    identical3(mpfr("0b101", base= 2),
               mpfr(  "101", base= 2), mpfr(5, precBits = 3))
   ,
    identical3(mpfr("0xabc", base=16),
               mpfr(  "abc", base=16), mpfr(2748, base=16, precBits = 12))
)


###----- _2_ Debugging, changing MPFR defaults, .. -----------------------------
##  NB: Currently mostly  *not* documented, not even .mpfr.erange()

stopifnot(Rmpfr:::.mpfr.debug() == 0,# the default level
          ## Debugging level 1:
          Rmpfr:::.mpfr.debug(1) == 0)# the previous level
r <- mpfr(7, 100)^-1000
r
## still prints as default

## where as this does print info: -- notably the very large values [3..6]:
.eranges <- function() sapply(names(Rmpfr:::.erange.codes), .mpfr.erange)
## now returning *double* - which loses some precision:
formatC(.eranges(), format="fg")

.mpfr.minPrec()
.mpfr.maxPrec()# debug printing shows the long integer (on 64 bit)

## Now, level 2 :
stopifnot(Rmpfr:::.mpfr.debug(2) == 1)
r
## with quite a bit of output

r2 <- r^100
r2
L <- r^-100000

## quite platform dependent {valgrind ==> bug? even in mpfr/gmp/.. ?}
str(.mpfr2list(x4))
x4 ## "similar info" as .mpfr2list(.)

## Increase maximal exponent:

tools:::assertWarning(
    .mpfr.erange.set("Emax", 5e18)) # too large {FIXME why only wann ??}
.mpfr.erange("Emax") # is unchanged
if(4e18 < .mpfr.erange("max.emax")) {
    .mpfr.erange.set("Emax", 4e18) # now ok:
    stopifnot(.mpfr.erange("Emax") == 4e18)
}


## revert to no debugging:
stopifnot(Rmpfr:::.mpfr.debug(0) == 2)
.mpfr.maxPrec()

L / (r2^-1000)# (could be more accurate?)

stopifnot(
    all.equal(L, r2^-1000, tol= 1e-27), # why not more accurate?
    all.equal(log(L), -100000 * (-1000) * log(7),
              tol = 1e-15)
)

## Now, our experimental "transport vehicle":
stopifnot(length(rv <- c(r, r2, L)) == 3)

str(mpfrXport(rv))
str(mpfrXport(mpfr(2, 64)^(-3:3)))
str(mpfrXport(Const("pi")* 2^(-3:3)))

## and a very large one
mil <- mpfr(1025, 111)
str(mm <- mpfrXport(xx <- mil^(2^25)))
stopifnot(all.equal(log2(xx) * 2^-25, log2(mil), tol=1e-15))

