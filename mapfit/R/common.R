## Common routines

diag.padding.zero <- function(A) {
  x <- diag(A)
  diag(A)[x == 0] <- NA
  A@x[is.na(A@x)] <- 0
  A
}

