outer.function <- function
### this is a high-level package function that will be exported for
### all of your users.
(a,
### first arg
 b
### second arg
 ) {
  foo <- a
  inner.function <- function
  ## This is a internal function only used inside outer.function, thus
  ## package users will not be able to access it. However we should
  ## still be able to document it like this.
  (x
   ## Some argument of the inner function.
   ) {
    y <- x
    ## Return value of the inner function.
  }
  bar <- b
### Return value of outer.function.
}
