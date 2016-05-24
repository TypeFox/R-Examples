.parseRelation <- function(x, y, rel) {
  return(switch(rel,
                "<="= x <= y,
                "<"= x < y,
                ">="= x >= y,
                ">"= x > y))
}

interval <- function(lower, upper, left=c(">=", ">"), right=c("<=", "<")) {
  left  <- match.arg(left)
  right <- match.arg(right)
  value <- list(lower=lower, upper=upper, left=left, right=right)
  class(value) <- "interval"
  return(value)
}

liesWithin <- function(x, int) {
  if (!check(int)) {
    warning("Not a valid interval!")
    return(rep(FALSE, length(x)))
  }
  return(sapply(x, function(y) {
    .parseRelation(y, int$lower, int$left) && .parseRelation(y, int$upper, int$right)
  }))
}

print.interval <- function(x, ...){
  leftBrace <- ifelse(x$left == ">", "(", "[")
  rightBrace <- ifelse(x$right == "<", ")", "]")
  interval <- paste(leftBrace, x$lower, ", ", x$upper, rightBrace, sep="")
  cat(interval, "\n")
  invisible(x)
}


check.interval <- function(object, ...) {
  isInterval <- TRUE
  leftOps  <- c(">", ">=")
  rightOps <- c("<", "<=")
  standardNames <- c("lower", "upper", "left", "right")
  haveSameElements <- identical(names(object), standardNames)
  if (!haveSameElements) {
    isInterval <- FALSE
    warning("The fields of the object do not equate the fields of an ",
            "'interval' object!", call.=FALSE)
    return(isInterval)
  }
  if (!object$left %in% leftOps) {
    isInterval <- FALSE
    warning("Invalid left comparison operator!", call.=FALSE)
    return(isInterval)
  }
  if (!object$right %in% rightOps) {
    isInterval <- FALSE
    warning("Invalid right comparison operator!", call.=FALSE)
    return(isInterval)
  }
  if (object$lower > object$upper) {
    isInterval <- FALSE
    warning("Empty interval!")
    return(isInterval)
  }
  return(isInterval)
}


