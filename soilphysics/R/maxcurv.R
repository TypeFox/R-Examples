maxcurv <- 
function (x.range, fun, 
	method = c("general", "pd", "LRP", "spline"), 
	x0ini = NULL, 
	graph = TRUE, ...) 
{
    stopifnot(is.atomic(x.range))
    if (!is.numeric(x.range)) 
        stop("'x.range' must be a numeric vector!")
    if (length(x.range) != 2) 
        stop("'x.range' must be a vector of length two!")
    if (diff(x.range) < 0) 
        stop("please, reorder 'x.range'.")
    if (!inherits(fun, "function")) 
        stop("'fun' must be a 'function' of x!")
    method <- match.arg(method)

    # first derivative
    dfun <- deriv3(fun2form(fun), "x", func = TRUE)
    if (attr(dfun(x.range[1]), "gradient") == attr(dfun(x.range[2]), 
        "gradient")) 
        stop("'fun' should not be a linar function of x!")

    # simulated points
    x <- seq(x.range[1], x.range[2], length.out = 5000)
    y <- fun(x)

    # ------------------------------------------------------
    # method: general curvature function (k) 
    if (method == "general") {
       # gradient and hessian
       gr <- attr(dfun(x), "gradient")
       he <- attr(dfun(x), "hessian")[,, "x"]

       # curvature function
       k <- abs(he) / (1 + gr^2)^(3/2)
       mcp <- x[which.max(k)]
    } else if (method == "pd") {
    # -------------------------------------------------------
    # method pd: perpendicular distances
       # parallel line
       b <- lm(range(fun(x.range)) ~ x.range)$coef
       if (fun(x.range[1]) > fun(x.range[2])) {
          si <- -1
       } else {
          si <- 1
       }
       b1 <- si * b[2]
       b0 <- mean(fun(x.range)) - b1 * mean(x.range)
       ang <- atan(b1)

       # perpendicular lines
       a1 <- -1/b1
       a0 <- y - a1*x
       xjs <- (b0 - a0) / (a1 - b1)
       yjs <- b0 + b1*xjs

       # perpendicular distances
       k <- sqrt( (x - xjs)^2 + (b0 + b1*xjs - y)^2 )
       mcp <- x[which.max(k)]
    } else if (method == "LRP") {
    # ---------------------------------------------------
    # method: LRP
       if (is.null(x0ini)) 
          stop("please, inform 'x0ini', a initial value for x0")
       ini <- coef(lm(y ~ x))
       fit <- try( nls(y ~ flrp(x, a0, a1, x0), 
             start = list(a0 = ini[1], a1 = ini[2], x0 = x0ini)), 
          silent = TRUE)
       if (class(fit) == "try-error") {
          fit <- try( nls(y ~ flrp(x, a0, a1, x0, left = FALSE), 
                start = list(a0 = ini[1], a1 = ini[2], x0 = x0ini)), 
             silent = TRUE)
       }
       if (class(fit) == "try-error") {
          stop("LRP could not get convergence!")
       } else {
          mcp <- coef(fit)[3]
          k <- rep(0, length(x))
       }
    } else {
    # ---------------------------------------------------
    # method: piecewise linear spline
       if (is.null(x0ini)) 
          stop("please, inform 'x0ini', a initial value for x0")
       lini <- coef(lm(y[x < x0ini] ~ x[x < x0ini]))
       rini <- coef(lm(y[x > x0ini] ~ x[x > x0ini]))
       fit <- try( nls(y ~ fs(x, a0, a1, b1, x0), 
             start = list(a0 = lini[1], a1 = lini[2], 
                b1 = rini[2], x0 = x0ini)), silent = TRUE)
       if (class(fit) == "try-error") {
          stop("spline could not get convergence!")
       } else {
          mcp <- coef(fit)[4]
          k <- rep(0, length(x))
       }
    }

    # graph
    if (graph) {
        curve(fun, from = x.range[1], to = x.range[2], ...)
        lines(x = c(mcp, mcp, -9e9), 
           y = c(-9e9, fun(mcp), fun(mcp)), lty = 3)
        if (method == "pd") {
           abline(b0, b1, lty = 3)
           lines(x = c(mcp, xjs[which.max(k)]),
              y = c(fun(mcp), b0 + b1*xjs[which.max(k)]), 
              lty = 3)
       } else if (method == "LRP" || method == "spline") {
          lines(x, predict(fit), col = 4, lty = 2)
       }
       devAskNewPage(ask = TRUE)
       plot(x, k, type = "l", ...)
       lines(x = c(mcp, mcp, -9e9), 
          y = c(-9e9, max(k), max(k)), lty = 3)
    }

    # output
    out <- list(fun = fun, 
       x0 = mcp, y0 = fun(mcp),
       method = method)
    class(out) <- "maxcurv"
    return(out)
}

# -------------------------------------------
# print method
print.maxcurv <- function (x, ...) 
{
    cat("\n          Maximum curvature point \n",
        "\nf(x) =", deparse(x$fun)[2],
        "\ncritical x: ", x$x0,
        "\ncritical y: ", x$y0)
    if (x$method == "general") {
       m <- "general curvature"
    } else if (x$method == "pd") {
       m <- "perpendicular distances"
    } else if (x$method == "LRP") {
       m <- "Linear Response Plateau"
    } else {
       m <- "piecewise linear spline"
    }
    cat("\nmethod:", m, "\n")
    invisible(x)
}

# -------------------------------------------
# LRP function
flrp <- function(x, a0, a1, x0, left = TRUE)
{
   y0 <- a0 + a1*x0
   if (left) {
      out <- ifelse(x <= x0, a0 + a1*x, y0)
   } else {
      out <- ifelse(x >= x0, a0 + a1*x, y0)
   }
   return(out)
}

# -------------------------------------------
# linear spline function
fs <- function(x, a0, a1, b1, x0)
{
  # left line: y = a0 + a1*x
  y0 <- a0 + a1*x0
  b0 <- y0 - (y0 - a0)*b1/a1
  ifelse(x <= x0, a0 + a1*x, b0 + b1*x)
}
