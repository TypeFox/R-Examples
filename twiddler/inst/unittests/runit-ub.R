##
## runit-ub.R - test for unboundVariables()
##

test.ub_01 <- function() {
  checkEquals(twiddler:::unboundVariables(quote(a + b)),
              list(quote(a), quote(b)))
}

test.ub_02 <- function() {
  a <- 1
  checkEquals(twiddler:::unboundVariables(quote(a + b)),
              list(quote(b)))
}

test.ub_03 <- function() {
  b <- 1
  checkEquals(twiddler:::unboundVariables(quote(a + b)),
              list(quote(a)))
}

test.ub_04 <- function() {
  a <- 1
  b <- 1
  checkEquals(twiddler:::unboundVariables(quote(a + b)),
              list())
}

## Make sure variables bound to functions (here c) are ignored and still returned as unbound
test.ub_04 <- function() {
  a <- 1
  b <- 1
  checkEquals(twiddler:::unboundVariables(quote(a + b + c)),
              list(quote(c)))
}

test.ubrec_01 <- function() {
  local({
    checkEquals(twiddler:::unboundVariables(quote(a + b)),
              list(quote(a), quote(b)))
  })
}

test.ubrec_02 <- function() {
  local({
    a <- 1
    checkEquals(twiddler:::unboundVariables(quote(a + b)),
                list(quote(b)))
  })
}

test.ubrec_03 <- function() {
  b <- 1
  local({
    checkEquals(twiddler:::unboundVariables(quote(a + b)),
                list(quote(a)))
  })
}

test.ubrec_04 <- function() {
  a <- 1
  local({
    b <- 1
    checkEquals(twiddler:::unboundVariables(quote(a + b)),
                list())
  })
}

## Trigger recursion using useless parentheses
test.ubrec_05 <- function() {
  local({
    b <- 1
    checkEquals(twiddler:::unboundVariables(quote(((((a)+0)+0)+0) + b)),
                list(quote(a)))
  })
}

test.ubfun_01 <- function() {
  checkEquals(twiddler:::unboundVariables(quote(f(a + b))),
              list(quote(a), quote(b)))  
}

test.ubfun_02 <- function() {
  b <- 1
  checkEquals(twiddler:::unboundVariables(quote(f(a + b))),
              list(quote(a)))  
}

test.ubfun_03 <- function() {
  a <- 1
  b <- 1
  checkEquals(twiddler:::unboundVariables(quote(f(a + b))),
              list())
}

## Ignore list like objects which are index into:
test.ubdf_01 <- function() {
  checkEquals(twiddler:::unboundVariables(quote(df$a + b)),
              list(quote(b)))
}

test.ubdf_02 <- function() {
  checkEquals(twiddler:::unboundVariables(quote(df[a] + b)),
              list(quote(a), quote(b)))
}

test.ubdf_03 <- function() {
  checkEquals(twiddler:::unboundVariables(quote(df[[a]] + b)),
              list(quote(a), quote(b)))
}
