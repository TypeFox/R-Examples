
## ************************************************************************
## Functions in this file are copied or modified from those in
## other R packages, as stated where appropriate.
## All authorship credit goes to the original author(s).
## 
## Note: A dot is prefixed to the original function names for internal usage.
## 
## 
## (c) Xiaobei Zhao
## 
## Wed Aug 06 16:26:27 EDT 2014 -0400 (Week 31)
## 
## 
## Reference: 
## 
## 
## ************************************************************************




## ------------------------------------------------------------------------
## functions from R package stringr
## ------------------------------------------------------------------------

.compact <- function(l) Filter(Negate(is.null), l)

.is.fixed <- function(string) {
  fixed <- attr(string, "fixed")
  if (is.null(fixed)) FALSE else fixed
}

.case.ignored <- function(string) {
  ignore.case <- attr(string, "ignore.case")
  if (is.null(ignore.case)) FALSE else ignore.case
}

.is.perl <- function(string) {
  perl <- attr(string, "perl")
  if (is.null(perl)) FALSE else perl
}

.recyclable <- function(...) {
  lengths <- vapply(list(...), length, integer(1))

  lengths <- lengths[lengths != 0]
  if (length(lengths) == 0) return(TRUE)

  all(max(lengths) %% lengths == 0)
}

.re_call <- function(f, string, pattern, replacement = NULL) {
  args <- list(pattern, replacement, string,
    fixed = .is.fixed(pattern), ignore.case = .case.ignored(pattern),
    perl = .is.perl(pattern))

  if (!("perl" %in% names(formals(f)))) {
    if (args$perl) message("Perl regexps not supported by ", f)
    args$perl <- NULL
  }

  do.call(f, .compact(args))
}


.check_string <- function(string) {
  if (!is.atomic(string))
    stop("String must be an atomic vector", call. = FALSE)

  if (!is.character(string))
    string <- as.character(string)

  string
}

.check_pattern <- function(pattern, string, replacement = NULL) {
  if (!is.character(pattern))
    stop("Pattern must be a character vector", call. = FALSE)

  if (!.recyclable(string, pattern, replacement)) {
    stop("Lengths of string and pattern not compatible")
  }

  pattern
}

.match_to_matrix <- function(match, global = FALSE) {
  if (global && length(match) == 1 && (is.na(match) || match == -1)) {
    null <- matrix(0, nrow = 0, ncol = 2)
    colnames(null) <- c("start", "end")

    return(null)
  }

  start <- as.vector(match)
  start[start == -1] <- NA
  end <- start + attr(match, "match.length") - 1L

  cbind(start = start, end = end)
}

.str_sub <- function(string, start = 1L, end = -1L) {
  if (length(string) == 0L || length(start) == 0L || length(end) == 0L) {
    return(vector("character", 0L))
  }

  string <- .check_string(string)

  n <- max(length(string), length(start), length(end))
  string <- rep(string, length = n)
  start <- rep(start, length = n)
  end <- rep(end, length = n)

  # Convert negative values into actual positions
  len <- .str_length(string)

  neg_start <- !is.na(start) & start < 0L
  start[neg_start] <- start[neg_start] + len[neg_start] + 1L

  neg_end <- !is.na(end) & end < 0L
  end[neg_end] <- end[neg_end] + len[neg_end] + 1L

  substring(string, start, end)
}


.str_length <- function(string) {
  string <- .check_string(string)

  nc <- nchar(string, allowNA = TRUE)
  is.na(nc) <- is.na(string)
  nc
}

.str_match <- function(string, pattern) {
  string <- .check_string(string)
  pattern <- .check_pattern(pattern, string)

  if (length(string) == 0) return(character())

  matcher <- .re_call("regexec", string, pattern)
  matches <- regmatches(string, matcher)

  # Figure out how many groups there are and coerce into a matrix with
  # nmatches + 1 columns
  tmp <- .str_replace_all(pattern, "\\\\\\(", "")
  n <- .str_length(.str_replace_all(tmp, "[^(]", "")) + 1

  len <- vapply(matches, length, integer(1))
  matches[len == 0] <- rep(list(rep(NA_character_, n)), sum(len == 0))

  do.call("rbind", matches)
}

.re_mapply <- function(f, string, pattern, replacement = NULL) {
  args <- list(
    FUN = f, SIMPLIFY = FALSE, USE.NAMES = FALSE,
    pattern, replacement, string,
    MoreArgs = list(
      fixed = .is.fixed(pattern),
      ignore.case = .case.ignored(pattern))
    )
  do.call("mapply", .compact(args))
}

.str_replace_all <- function(string, pattern, replacement) {
  string <- .check_string(string)
  pattern <- .check_pattern(pattern, string, replacement)

  if (length(pattern) == 1 && length(replacement) == 1) {
    .re_call("gsub", string, pattern, replacement)
  } else {
    unlist(.re_mapply("gsub", string, pattern, replacement))
  }
}


.str_locate_all <- function(string, pattern) {
  string <- .check_string(string)
  pattern <- .check_pattern(pattern, string)

  if (length(pattern) == 1) {
    matches <- .re_call("gregexpr", string, pattern)
  } else {
    matches <- unlist(.re_mapply("gregexpr", string, pattern),
      recursive = FALSE)
  }
  lapply(matches, .match_to_matrix, global = TRUE)
}


.str_extract_all <- function(string, pattern) {
  string <- .check_string(string)
  pattern <- .check_pattern(pattern, string)

  positions <- .str_locate_all(string, pattern)
  lapply(seq_along(string), function(i) {
    position <- positions[[i]]
    .str_sub(string[i], position[, "start"], position[, "end"])
  })
}


.str_match_all <- function(string, pattern) {
  matches <- .str_extract_all(string, pattern)

  lapply(matches, function(match) {
    .str_match(match, pattern)
  })
}
