#### Experiment with spline code --- want to use
library(splines) ## for the splines
suppressMessages(library(cobs))

options(digits = 9)

### -- 1) -- Look at `splines' pkg code :
data(women)
yw <- women$weight
xh <- women$height# == 58:72; too trivial; modify a bit:
ii <- c(2,5,9,11); xh[ii] <- xh[ii]+ 0.2
ii <- c(3:4,8,13); xh[ii] <- xh[ii]+ 0.25
(bIspl <- interpSpline( xh, yw, bSpl = TRUE))
print.default(bIspl[1:3])
(pIspl <- interpSpline( xh, yw, bSpl = FALSE))
p2Ispl <- polySpline(bIspl)
all.equal(pIspl, p2Ispl, tol = 1e-15)# TRUE
##--> could use polySpline() at end of interpSpline(.)


### -- 2) --- substituting our .splBasis()  by splines package splineDesign() --
### ========
### ---> done (in principle; not yet implemented!),  Feb.2002

str(cobs:::.splBasis(4, bIspl$knots, length(bIspl$coef) + 6, x = .5 + 57:72))# outside!
xo <- 0.5 + 59:70 # should work up to ord = 5

## ord <- 4 # cubic splines
## ord <- 3 # quadratic splines
for(ord in 5:1) {
    cat("\n\nord = ",ord,"\n========\n")
    print(spB <- cobs:::.splBasis(ord, bIspl$knots,
                           length(bIspl$coef) + ord + 2, x = xo))
    ## Gives error for ord = 5:4 --- data must be INSIDE :
    try(       splineDesign(bIspl$knots, x = 0.5 + 57:72, ord = ord))
    str(spD <- splineDesign(bIspl$knots, x = xo,          ord = ord))

    ## splineDesign() contains:
    tmp <- .Call("spline_basis", bIspl$knots, ord=ord, x= xo,
                 derivs = integer(length(xo)), PACKAGE = "splines")
    print(offs.tmp <- attr(tmp, "Offsets"))
    attr(tmp, "Offsets") <- NULL
    print(all.equal(tmp, spB$design, tol = 4e-16)) # TRUE
    print(all(spB$offsets - offs.tmp == ord - 1)) # TRUE
}

## This (comments dropped) checks the same but gives no output (iff OK)
for(ord in 5:1) {
    spB <- cobs:::.splBasis(ord, bIspl$knots, length(bIspl$coef) + ord + 2, x = xo)
    tmp <- .Call("spline_basis", bIspl$knots, ord=ord, x= xo,
                 derivs = integer(length(xo)), PACKAGE = "splines")
    offs.tmp <- attr(tmp, "Offsets")
    attr(tmp, "Offsets") <- NULL
    stopifnot(all.equal(tmp, spB$design, tol = 4e-16),
              all(spB$offsets - offs.tmp == ord - 1))
}

### -- 3) --- substituting our cobs:::.splValue() by splines package predict.bSpline
### ========
### ----------- STILL TODO !! ------------

### ~/R/D/r-devel/R/src/library/splines/R/splineClasses.R  has
##    predict.polySpline
##    predict.bSpline    <- function(object, x, nseg = 50, deriv = 0, ...)
##    predict.nbSpline
##    predict.pbSpline         ((all with the same argument list))
##    predict.npolySpline
##    predict.ppolySpline

## Using  bIspl ,  the interpolating B-spline from above
predict(bIspl, xo)
## $x
##  [1] 59.5 60.5 61.5 62.5 63.5 64.5 65.5 66.5 67.5 68.5 69.5 70.5
##
## $y
##  [1] 117.761938 120.761884 123.743134 127.112180 130.660049 133.066681
##  [7] 135.940123 140.232337 143.472651 147.532222 151.495247 155.519340
##
## attr(,"class")
## [1] "xyVector"

cobs:::.splValue(deg = 3, knots = bIspl$knots, coef = bIspl$coef, xo = xo)
## hmm, not the same as $ y above ...
