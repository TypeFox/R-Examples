##
##  r e g e x p . R
##

regexp <- function(s, pat, ignorecase = FALSE, once = FALSE, split = FALSE)
{
    # Match regular expression
    stopifnot(is.character(pat), is.character(s))
    if (length(pat) > 1) {
        warning("Only the first string in argument 'pat' is taken.")
        pat <- pat[1]
    }
    if (length(s) > 1) {
        warning("Only the first string in argument 's' is taken.")
        s <- s[1]
    }
    if (is.na(pat) || is.na(s))
        stop("In arguments 'pat' and 's' NA values not allowed.")

    if (once) {
        res <- regexpr(pat, s, ignore.case = ignorecase, perl = TRUE)
    } else {
        res <- gregexpr(pat, s, ignore.case = ignorecase, perl = TRUE)[[1]]
    }

    if (length(res) == 1 && res < 0)
        if (split)
            return(list(start = NULL, end = NULL, match = NULL, split = s))
        else
            return(list(start = NULL, end = NULL, match = NULL, split = NULL))

    # Return start, end, and match
    rstart <- res
    rend <- rstart + attr(res, "match.length") - 1
    attr(rstart, "match.length") <- attr(rend, "match.length") <-NULL
    rmatch <- substring(s, rstart, rend)
    if (split) {
        # rsplit <- strsplit(s, pat, perl = TRUE)[[1]]  # does not ignore case
        n <- nchar(s)
        rs <- c(0, rstart, n+1)
        re <- c(0, rend, n+1)
        rsplit <- c()
        for (i in 1:(length(rs)-1)) {
            if (rs[i+1] - re[i] > 1)
                rsplit <- c(rsplit, substr(s, re[i]+1, rs[i+1]-1))
        }
    } else {
        rsplit <- NULL
    }

    list(start = rstart, end = rend, match = rmatch, split = rsplit)
}

regexpi <- function(s, pat, once = FALSE, split = FALSE)
{
    regexp(s, pat, ignorecase = TRUE, once = once, split = split)
}

regexprep <- function(s, expr, repstr, ignorecase = FALSE, once = FALSE)
{
    # Replace string using regular expression
    if (! is.character(s))
        stop("Argument 's' must be a character vector.")
    if (!is.character(expr) || !is.character(repstr) ||
        length(expr) != 1   || length(repstr) != 1)
        stop("Arguments 'old' and 'new' must be simple character strings.")

    if (once) {
        sub(expr, repstr, s, ignore.case = ignorecase, perl = TRUE)
    } else {
        gsub(expr, repstr, s, ignore.case = ignorecase, perl = TRUE)
    }
}

refindall <- function(s, pat, over = 1, ignorecase = FALSE)
{
  if (!is.character(s) || !is.character(pat) ||
      length(s) != 1 || length(pat) != 1)
      stop("Arguments 's' and 'pat' must be single strings.")
  if (!is.numeric(over) || length(over) != 1 ||
      over < 1 || over != round(over))
      stop("Argument 'over' must be a positive integer.")

  pos <- c()  # positions of matches
  i <- 1; n <- nchar(s)
  found <- regexpr(pat, substr(s, i, n), ignore.case = ignorecase, perl=TRUE)
  while (found > 0) {
    pos <- c(pos, i + found - 1)
    i <- i + found + (over - 1)
    found <- regexpr(pat, substr(s, i, n), ignore.case = ignorecase, perl=TRUE)
  }
  return(pos)
}
