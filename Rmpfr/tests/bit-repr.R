stopifnot(suppressPackageStartupMessages(require("Rmpfr")))
## (checking that the 32 / 64 bit  GMP message does *not* show here)

### Try to look at the internal bit-representation of the limbs

.limbs <- function(x) {
    stopifnot(is(x, "mpfr"))
    lapply(x@.Data, slot, "d") # not sapply() each can have different prec. & #{limbs}
}
.expo <- function(x) {
    stopifnot(is(x, "mpfr"))
    sapply(x@.Data, slot, "exp")
}

Bits <- function(x) {
    L <- .limbs(x)# list(length n) each of "k(prec)" 32-bit ints
    ## NB:  mpfr(2, .) and all mpfr(2^k, .) also have a 'd' ending in NA integer!
    ##     [reason: after all, R's NA_integer_ is INT_MAX+1 = 2^31 ]
    ## and  the  mpfr(c(NA,NaN, Inf, -Inf), .)   have *no* NA in 'd' (but all in 'exp'!
    ## see .mpfr2list() example below

    hasNA <- any(iNA <- sapply(lapply(L, is.na), any)) # iNA: TRUE if there's an NA
    ## need to catch them later
    CC <- function(ch) paste(ch, collapse="")
    hex <- sapply(L, function(.) CC(sprintf("%x", rev(.))))
    if(hasNA) hex[iNA] <- NA_character_
    hex <- strsplit(hex, NULL)

    db <- t(expand.grid(0:1,0:1,0:1,0:1, KEEP.OUT.ATTRS=FALSE)[,4:1])
    storage.mode(db) <- "character" # "0" or "1"
    dimnames(db) <- list(NULL, c(paste(0:9), letters[1:6]))
    ## db is  4 x 16  matrix  with col.names "0" "1" .. "9" "a" "b" ... "f"

    ex <- .expo(x)
    if(is.matrix(ex)) {
        ## 64-bit case: exponent is long == two ints
        ## -----------  the 2nd int is in {0, -1, NA} (NA : for 0)
        ex2 <- ex[2,]
        ex <- ex[1,]
    }
    pat <- paste("(", sapply(pmax(0, ex),
                             function(n) CC(rep.int(".", n))),
                 ")0+$", sep="")
    ## pat <- ifelse(iNA, NA_character_, pat)

    getbits <- function(ch) CC(as.vector(db[,ch]))

    ## drop trailing zeros (from r[[i]], via pat[i]), keeping correct number:
    drop0.r <- function(i) sub(pat[i], "\\1", r[[i]])

    if(hasNA) {
        r <- as.list(iNA)
        r[!iNA] <- lapply(hex[!iNA], getbits)
        r[!iNA] <- lapply(which(!iNA), drop0.r)
        ## FIXME this is wrong -- really have powers of 2, and want their (easy) bits :
        r[iNA ] <- NA_character_
        unlist(r)
    }
    else {
	r <- lapply(hex, getbits)
        sapply(seq_along(r), drop0.r)
    }

}

x <- mpfr(r <- c(NA,NaN, Inf, -Inf), 64)
stopifnot(identical(asNumeric(x), # mpfr has no NA, just NaN's:
                    c(NaN,NaN, Inf, -Inf)))

if(FALSE) # platform dependent:
    ## The "non-finite" mpfr value internals (in 64-bit: 'exp' has NA):
    str(.mpfr2list(x))


## bug in Bits(): all (exact) powers of 2  will show as NA:

x <- mpfr(c(3:5,11:16, 59, 125:128, 1024:1025), 64)
x
data.frame(x= as.numeric(x), I(Bits(x)))

x <- mpfr(c(-20:-1, 1:30),64)# w/o  0  - as its mantissa is "random" (in 64-bit)
data.frame(x= as.numeric(x), I(Bits(x)))
b0 <- Bits(mpfr(0, 64)) # not printing it here -- they are "random" for this special case!

(half <- mpfr(0.5, 64)*(1 + mpfr(2, 64)^-16 * (-3:3)))
Bits(half)

## pi, in varying number of bits :
p. <- round(pi* 2^c(10,16,5*(4:8)))
dput(p.)#-> the definition of p :
p <- mpfr(c(3217, 205887, 3294199, 105414357,
            3373259426, 107944301636, 3454217652358), 64)
stopifnot(all.equal(p., p, tolerance = 1e-15))
## all the mantissas are those of pi, rounded differently:
Bits(c(p, Const("pi", 64)))


cat('Time elapsed: ', proc.time(),'\n') # "stats"

if(!interactive()) warnings()
