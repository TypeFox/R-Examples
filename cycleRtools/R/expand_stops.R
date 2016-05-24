## Simple example:
# d <- data.frame(timer.s = c(1:50, 101:150), x = rnorm(100))
# d <- expand_stops(d); attributes(d)
expand_stops <- function(data, deltat = NULL, tcol = "timer.s",
                         keep_attr = FALSE) {
  tdat  <- data[, tcol]
  diffs <- table(Diff(tdat))
  diffs <- rbind(count = unname(diffs), val = as.numeric(names(diffs)))

  if (is.null(deltat)) {
    # deltat should represent > 75% of sample time differences.
    if ((diffs["count", which.max(diffs["count", ])] / length(tdat)) < 0.75)
      stop(bad_deltat_msg(), call. = FALSE)
    deltat <- diffs["val", which.max(diffs["count", ])]
  }
  if (deltat != 1 && deltat != 0.5)
    warning(odd_deltat_msg(), call. = FALSE)

  # "stop" = difference between samples > 10 sec.
  stops <- diffs["val", diffs["val", ] > 10]
  stopi <- which(c(0, Diff(tdat)) %in% stops)

  fill  <- unlist(lapply(stopi, function(i)
    seq(from = tdat[i - 1], to = tdat[i], by = deltat)))
  fill  <- setdiff(fill, data[, tcol])

  # No stops to expand.
  if (is.null(fill)) {
    attr(data, "wo_expand") <- seq_along(tdat)
    attr(data, "deltat")    <- unname(deltat)
    return(data)
  }

  NAs   <- rep_len(NA, length(fill))
  empty <- lapply(colnames(data), function(x) NAs)
  empty <- setNames(empty, colnames(data))
  empty[[tcol]] <- fill
  empty <- as.data.frame(empty)

  out <- merge(data, empty, by = tcol, all = TRUE)
  out <- setNames(out[, c(1, grep("\\.x$", colnames(out)))], colnames(data))

  if (keep_attr)
    attributes(out) <- attributes(data)

  # Metadata to be used elsewhere.
  attr(out, "new")       <- match(fill, out[, tcol])
  attr(out, "wo_expand") <- which(out[, tcol] %notin% fill)
  attr(out, "deltat")    <- unname(deltat)

  out
}

bad_deltat_msg <- function()
  paste("Ambiguous sampling frequency. Data should be sampled uniformly.",
        "Revise data and/or supply deltat argument.", sep = "\n")
odd_deltat_msg <- function()
  paste("Unusual deltat value produced; output may be spurious.",
        "Make sure data is sampled *regularly*.", sep = "\n")
