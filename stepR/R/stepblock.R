"stepblock" <-
function(value, leftEnd = c(1, rightEnd[-length(rightEnd)] + 1), rightEnd = 1:length(value), x0 = 0)
{
  ret <- data.frame(leftEnd = leftEnd, rightEnd = rightEnd, value = value)
  attr(ret, "x0") <- x0
  class(ret) <- c("stepblock", class(ret))
  ret
}

"[.stepblock" <- 
function (x, i, j, drop = if(missing(i)) TRUE else if(missing(j)) FALSE else length(j) == 1, ...) 
{
  ret <- as.data.frame(x)
  ret <- "[.data.frame"(ret, i, j, drop)
  if(!missing(i)) { # if no rows are selected, return a data.frame for convenience
    a <- attributes(x)
    if(!missing(i)) {
      i <- (1:nrow(x))[i]
      a$row.names <- a$row.names[i]
      ret$leftEnd <- x$leftEnd[c(0, i[-length(i)]) + 1]
      ret$leftIndex <- x$leftIndex[c(0, i[-length(i)]) + 1]
    }
    if(!missing(j)) a$names <- names(ret)
    attributes(ret) <- a
  }
  ret
}

"print.stepblock" <-
function(x, ...)
{
  cat("\n")
  cat("Step function containing", nrow(x), "blocks\n\n")
  cat("domain: (", attr(x, "x0"), ",", x$rightEnd[nrow(x)], "]\n")
  cat("range:  [", min(x$value), ",", max(x$value), "]\n")
  cat("\n")
}

"plot.stepblock" <-
function(x, type = "c", xlab = "x", ylab = "y", main = "Step function", sub = NULL, ...)
{
  xe <- c(attr(x, "x0"), x$rightEnd)
  xb <- c(x$leftEnd, x$rightEnd[nrow(x)])
  ys <- c(x$value[1], x$value)
  newtype <- "S"
  switch(tolower(type),
    c = {
      if(is.null(sub)) sub <- "jumps between blocks"
      xs <- ( xe + xb ) / 2
    },
    e = {
      if(is.null(sub)) sub <- "jumps at end of block"
      xs <- xe
    },
    b = {
      if(is.null(sub)) sub <- "jumps at beginning of block"
      xs <- xb
    },
    {
      if(is.null(sub)) sub <- ""
      newtype <- type
      xs <- x$x[x$rightEnd]
      ys <- x$value
    }
  )
  plot(xs, ys, type = newtype, xlab = xlab, ylab = ylab, main = main, sub = sub, ...)
  switch(EXPR = type,
    C = {
      if(is.null(match.call()$pche)) pche <- 1
      if(is.null(match.call()$pche)) pchb <- 1
      pointsbe <- T
    },
    E = {
      if(is.null(match.call()$pche)) pche <- 19
      if(is.null(match.call()$pche)) pchb <- 1
      pointsbe <- T
    },
    B = {
      if(is.null(match.call()$pche)) pche <- 1
      if(is.null(match.call()$pche)) pchb <- 19
      pointsbe <- T
    },
    pointsbe <- F
  )
  if(pointsbe) {
    points(xs[-1], x$value, pch = pche, ...)
    points(xs[-length(xs)], x$value, pch = pchb, ...)
  }
}

"lines.stepblock" <-
function(x, type = "c", ...)
{
  xe <- c(attr(x, "x0"), x$rightEnd)
  xb <- c(x$leftEnd, x$rightEnd[nrow(x)])
  ys <- c(x$value[1], x$value)
  newtype <- "S"
  switch(tolower(type),
    c = {
      xs <- ( xe + xb ) / 2
    },
    e = {
      xs <- xe
    },
    b = {
      xs <- xb
    },
    {
      newtype <- type
      xs <- x$x[x$rightEnd]
      ys <- x$value
    }
  )
  lines(xs, ys, type = newtype, ...)
  if(!"pche" %in% names(match.call())) {
    if(type == "E") pche <- 19 else pche <- 1
  } else pche <- eval(match.call()[["pche"]])
  if(!"pchb" %in% names(match.call())) {
    if(type == "B") pchb <- 19 else pchb <- 1
  } else pchb <- eval(match.call()[["pche"]])
  if(type %in% c("C", "E", "B")) {
    points(xs[-1], x$value, pch = pche, ...)
    points(xs[-length(xs)], x$value, pch = pchb, ...)
  }
}
