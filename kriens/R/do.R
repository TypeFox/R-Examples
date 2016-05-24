do <- function(...) {
  args <- list(...)
  # prevents lazy evalutation
  force(args)
  forget(path(args))
}
