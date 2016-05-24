path <- function(fs) {
  # prevents lazy evaluation
  force(fs)
  if(length(fs) == 0 || is.null(fs))
    stop("the provided list of function is not valid (either empty or null)")
  Reduce(compose, fs)
}
