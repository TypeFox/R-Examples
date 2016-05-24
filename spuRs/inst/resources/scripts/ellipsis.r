test4 <- function(x, ...) {
  return(sd(x, ...))
}

test4(1:3)

test4(c(1,NA,3))

test4(c(1,NA,3), na.rm = TRUE)

test5 <- function(x, na.rm = FALSE, ...) {
  return(sd(x, ...))
}

test5(c(1,NA,3))

test5(c(1,NA,3), na.rm = TRUE)

test6 <- function(x, na.rm = FALSE, ...) {
  return(sd(x, na.rm = na.rm, ...))
}

test6(c(1,NA,3))

test6(c(1,NA,3), na.rm = TRUE)



