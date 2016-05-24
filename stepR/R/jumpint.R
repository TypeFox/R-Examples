"jumpint" <-
function (sb, ...) UseMethod("jumpint")

"jumpint.stepfit" <-
function (sb, ...)
{
  ci.names <- c("leftEndLeftBound", "leftEndRightBound", "rightEndLeftBound", "rightEndRightBound",
    "leftIndexLeftBound", "leftIndexRightBound", "rightIndexLeftBound", "rightIndexRightBound")
  if(!all(ci.names %in% names(sb))) stop("sb does not contain confidence intervals")
  ci <- as.data.frame(sb)[,ci.names]
  attr(ci, "x0") <- attr(sb, "x0")
  class(ci) <- c("jumpint", class(ci))
  ci
}

"points.jumpint" <-
function(x, pch.left = NA, pch.right = NA, y.left = NA, y.right = NA, xpd = NA, ...)
{
  if(is.na(pch.left)) pch.left <- "("
  if(is.na(pch.right)) pch.right <- "]"
  if(is.na(y.left)) y.left <- par()$usr[3]
  if(is.na(y.right)) y.right <- par()$usr[3]
  points(x$rightEndLeftBound[-nrow(x)], if(length(y.left) == 1) rep(y.left, nrow(x) - 1) else y.left, pch = pch.left, xpd = xpd, ...)
  points(x$leftEndRightBound[-1], if(length(y.right) == 1) rep(y.right, nrow(x) - 1) else y.right, pch = pch.right, xpd = xpd, ...)
}

"confband" <-
function (sb, ...) UseMethod("confband")

"confband.stepfit" <-
function (sb, ...)
{
  if(is.null(attr(sb, "confband")))  stop("sb does not contain confidence bands")
  cb <- attr(sb, "confband")
  for(a in c("family", "param")) attr(cb, a) <- attr(sb, a)
  cb
}

"lines.confband" <-
function(x, dataspace = TRUE, ...)
{
  scale <- if(attr(x, "family") == "binomial" && dataspace) attr(x, "param") else 1
  lines(c(x$x, x$x[nrow(x)]), scale * c(x$lower[1], pmin(x$lower, c(x$lower[-1], Inf))), type = "S", ...)
  lines(c(x$x[1], x$x), scale * c(pmax(x$upper, c(-Inf, x$upper[-nrow(x)])), x$upper[nrow(x)]), type = "s", ...)
}
