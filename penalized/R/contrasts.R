#######################################
# An alternative contrasts function for unordered factors
# Ensures symmetric treatment of all levels of a factor
#######################################
contr.none <- function(n, contrasts) {
  if (length(n) == 1)
    contr.treatment(n, contrasts = n<=2)
  else
    contr.treatment(n, contrasts = length(unique(n))<=2)
}

#######################################
# An alternative contrasts function for ordered factors
# Ensures use of a difference penalty for such factors
#######################################
contr.diff <- function (n, contrasts = TRUE)
{
    if (is.numeric(n) && length(n) == 1) {
        if (n > 1)
            levs <- 1:n
        else stop("not enough degrees of freedom to define contrasts")
    }
    else {
        levs <- n
        n <- length(n)
    }
    contr <- array(0, c(n, n), list(levs, paste(">=", levs, sep="")))
    contr[outer(1:n,1:n, ">=")] <- 1
    if (n < 2)
      stop(gettextf("contrasts not defined for %d degrees of freedom",
        n - 1), domain = NA)
    if (contrasts)
      contr <- contr[, -1, drop = FALSE]
    contr
}