
rvvapply <- function (FUN, n., ..., constantArgs = list()) {
  args <- list(...)
  n.sims <- getnsims()
  if (length(constantArgs) > 0 && any(sapply(constantArgs,
                                             is.rv))) {
    stop("constantArgs should contain only constant arguments")
  }
  max_vector_length <- max(sapply(args, length))
  .resize_sims <- function(x, ...) {
    s <- t(sims(as.rv(x)))
    resizeSims(s, ...)
  }
  args <- lapply(args, .resize_sims, vector.length = max_vector_length,
                 n.sims = n.sims)
  n.name <- attr(n., "n.name")
  if (is.null(n.name)) {
    n.name <- "n"
  }
  random_n <- (!missing(n.) && is.random(n.))
  if (!missing(n.)) {
    if (length(n.) > 1) {
      stop("length(n.)>1? Random dimensions not supported")
    }
    n. <- n.[1]
    n.max <- rvmax(n.)
    args[[n.name]] <- (n.max * max_vector_length * n.sims)
  }
  args <- c(args, constantArgs)
  FUN <- match.fun(FUN)
  S <- do.call(FUN, args = args)
  n.max <- (length(S) %/% (max_vector_length*n.sims))
  S <- array(S, c(max_vector_length, n.sims, n.max))
  S <- aperm(S, c(2, 1, 3))
  S <- t(matrix(S, nrow=n.sims))
  if (random_n) {
    # Here we must have variables in rows, simulations in columns
    n.max <- rvmax(n.)
    ix <- rep(1:n.max, each = max_vector_length)
    ns <- as.vector(sims(n.))
    not_observed <- sapply(ns, function(n) ix > n)
    S[not_observed] <- NA
  }
  x <- rvsims(t(S))
  n.scalars <- length(x)
  n.rows <- (n.scalars %/% max_vector_length)
  if (n.rows > 1) {
    if (max_vector_length == 1) {
      dim(x) <- NULL
    }
    else {
      dim(x) <- c(max_vector_length, n.rows)
    }
  }
  return(x)
}
