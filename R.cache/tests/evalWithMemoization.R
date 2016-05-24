library("R.cache")

for (kk in 1:5) {
  cat(sprintf("Iteration #%d:\n", kk))
  res <- evalWithMemoization({
    cat("Evaluating expression...")
    a <- 1
    b <- 2
    c <- 4
    cat("done\n")
    b
  })
  print(res)

  # Sanity checks
  stopifnot(a == 1 && b == 2 && c == 4)

  # Clean up
  rm(list=c("a", "b", "c"))
} # for (kk ...)


############################################################
# WARNING
############################################################
# If the expression being evaluated depends on
# "input" objects, then these must be be specified
# explicitly as "key" objects.
for (ii in 1:2) {
  for (kk in 1:3) {
    cat(sprintf("Iteration #%d:\n", kk))
    res <- evalWithMemoization({
      cat("Evaluating expression...")
      a <- kk
      cat("done\n")
      a
    }, key=list(kk=kk))
    print(res)

    # Sanity checks
    stopifnot(a == kk)

    # Clean up
    rm(list=c("a"))
  } # for (kk ...)
} # for (ii ...)
