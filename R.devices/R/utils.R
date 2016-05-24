nullfile <- function() {
  switch(.Platform$OS.type,
    windows="NUL",
    "/dev/null"
  )
}
