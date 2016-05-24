##
## cec2009 unit test
##

## Source to generate test point checks:
## 
## f <- function(n, x=c(0.2, 0.2, 0.2, 0.2)) {
##   fn <- get(n)
##   r <- fn(x)
##   message("checkEqualsNumeric(", n, "(",
##           "c(", paste(x, collapse=", "), "), ", 
##           "c(", paste(r, collapse=", "), ")))")
## }
##
## for (n in paste("UF", 1:10, sep=""))
##   f(n)



test.UF1 <- function() {
  checkEqualsNumeric(UF1(c(0.2, 0.2, 0.2, 0.2)), c(0.454091055737032, 1.72127910133303))
  checkEqualsNumeric(UF1(c(1.0, 0.0, 0.0)), c(1.0, 1.5))
  checkEqualsNumeric(UF1(c(1.0, 2.0, 2.0)), c(NaN, NaN))
  checkEqualsNumeric(UF1(c(-1.0, 0.0, 0.0)), c(NaN, NaN))
  checkEqualsNumeric(UF1(c(0.2, 0.2, "a")), c(NaN, NaN))
  checkException(UF1(c(1.0, 2.0)))
}

test.UF2 <- function() {
  checkEqualsNumeric(UF2(c(0.2, 0.2, 0.2, 0.2)), c(0.210335975924134, 0.654710703252405))  
  checkEqualsNumeric(UF2(c(1.0, 1.0, 1.0)), c(8.2199999999999989, 0.7449042731880116))
  checkEqualsNumeric(UF2(c(1.0, 2.0, 2.0)), c(NaN, NaN))
  checkEqualsNumeric(UF2(c(-1.0, 0.0, 0.0)), c(NaN, NaN))
  checkEqualsNumeric(UF2(c(0.2, 0.2, "a")), c(NaN, NaN))
  checkException(UF2(c(1.0, 2.0)))
}

test.UF3 <- function() {
  checkEqualsNumeric(UF3(c(0.2, 0.2, 0.2, 0.2)), c(7.1937361619286, 2.90716200450771))
  checkEqualsNumeric(UF3(c(1.0, 2.0, 2.0)), c(NaN, NaN))
  checkEqualsNumeric(UF3(c(-1.0, 0.0, 0.0)), c(NaN, NaN))
  checkEqualsNumeric(UF3(c(0.2, 0.2, "a")), c(NaN, NaN))
  checkException(UF3(c(1.0, 2.0)))
}

test.UF4 <- function() {
  checkEqualsNumeric(UF4(c(0.2, 0.2, 0.2, 0.2)), c(0.434509085633295, 1.20063945204195))
  checkEqualsNumeric(UF4(c(1.0, 2.0, 2.0)), c(NaN, NaN))
  checkEqualsNumeric(UF4(c(-1.0, 0.0, 0.0)), c(NaN, NaN))
  checkEqualsNumeric(UF4(c(0.2, 0.2, "a")), c(NaN, NaN))
  checkException(UF4(c(1.0, 2.0)))
}

test.UF5 <- function() {
  checkEqualsNumeric(UF5(c(0.2, 0.2, 0.2, 0.2)), c(3.17056357285719, 3.98342430089429))
}

test.UF6 <- function() {
  checkEqualsNumeric(UF6(c(0.2, 0.2, 0.2, 0.2)), c(1.88938265125039, 6.65748253841339))
}

test.UF7 <- function() {
  checkEqualsNumeric(UF7(c(0.2, 0.2, 0.2, 0.2)), c(0.978870719414727, 1.44371303315529))
}

test.UF8 <- function() {
  checkEqualsNumeric(UF8(c(2.0, 0.2, 0.2, 0.2, 0.2)), c(NaN, NaN, NaN))
  checkEqualsNumeric(UF8(c(0.2, 2.0, 0.2, 0.2, 0.2)), c(NaN, NaN, NaN))
  checkEqualsNumeric(UF8(c(0.2, 0.2, 3.0, 0.2, 0.2)), c(NaN, NaN, NaN))
  checkEqualsNumeric(UF8(c(0.2, 0.2, 0.2, 0.2, "a")), c(NaN, NaN, NaN))
  for (i in 1:4)
    checkException(UF8(rep(0.2, i)))
}

test.UF9 <- function() {
  checkEqualsNumeric(UF9(c(2.0, 0.2, 0.2, 0.2, 0.2)), c(NaN, NaN, NaN))
  checkEqualsNumeric(UF9(c(0.2, 2.0, 0.2, 0.2, 0.2)), c(NaN, NaN, NaN))
  checkEqualsNumeric(UF9(c(0.2, 0.2, 3.0, 0.2, 0.2)), c(NaN, NaN, NaN))
  checkEqualsNumeric(UF9(c(0.2, 0.2, 0.2, 0.2, "a")), c(NaN, NaN, NaN))
  for (i in 1:4)
    checkException(UF9(rep(0.2, i)))
}

test.UF10 <- function() {
  checkEqualsNumeric(UF10(c(2.0, 0.2, 0.2, 0.2, 0.2)), c(NaN, NaN, NaN))
  checkEqualsNumeric(UF10(c(0.2, 2.0, 0.2, 0.2, 0.2)), c(NaN, NaN, NaN))
  checkEqualsNumeric(UF10(c(0.2, 0.2, 3.0, 0.2, 0.2)), c(NaN, NaN, NaN))
  checkEqualsNumeric(UF10(c(0.2, 0.2, 0.2, 0.2, "a")), c(NaN, NaN, NaN))
  for (i in 1:4)
    checkException(UF10(rep(0.2, i)))
}
