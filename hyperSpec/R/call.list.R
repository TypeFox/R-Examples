###-----------------------------------------------------------------------------
###
### generate a list of function arguments for the calling function
###

.call.list <- function (x = NULL) {
  if (is.null (x))
    x <- sys.call (-1)

  if (length (x) < 3L)
    I (list ())
  else {
    x <- as.list (x [- (1 : 2)])
    I (x)
  }
}
