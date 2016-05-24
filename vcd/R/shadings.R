## convenience function for interfacing
## HCL colors as implemented in colorspace
hcl2hex <- function(h = 0, c = 35, l = 85, fixup = TRUE)
{
  colorspace::hex(polarLUV(l, c, h), fixup = fixup)
}

## shading-generating functions should take at least the arguments
##   observed, residuals, expected, df
## and return a function which takes a single argument (interpreted
## to be a vector of residuals).

shading_hsv <- function(observed, residuals = NULL, expected = NULL, df = NULL,
  h = c(2/3, 0), s = c(1, 0), v = c(1, 0.5),
  interpolate = c(2, 4), lty = 1, eps = NULL, line_col = "black",
  p.value = NULL, level = 0.95, ...)
{
  ## get h/s/v and lty
  my.h <- rep(h, length.out = 2)  ## positive and negative hue
  my.s <- rep(s, length.out = 2)  ## maximum and minimum saturation
  my.v <- rep(v, length.out = 2)  ## significant and non-significant value
  lty <- rep(lty, length.out = 2) ## positive and negative lty

  ## model fitting (if necessary)
  if(is.null(expected) && !is.null(residuals)) stop("residuals without expected values specified")
  if(!is.null(expected) && is.null(df) && is.null(p.value)) {
    warning("no default inference available without degrees of freedom")
    p.value <- NA
  }
  if(is.null(expected) && !is.null(observed)) {
    expected <- loglin(observed, 1:length(dim(observed)), fit = TRUE, print = FALSE)
    df <- expected$df
    expected <- expected$fit
  }
  if(is.null(residuals) && !is.null(observed)) residuals <- (observed - expected)/sqrt(expected)

  ## conduct significance test (if specified)
  if(is.null(p.value)) p.value <- function(observed, residuals, expected, df)
    pchisq(sum(as.vector(residuals)^2), df, lower.tail = FALSE)
  if(!is.function(p.value) && is.na(p.value)) {
    v <- my.v[1]
    p.value <- NULL
  } else {
    if(is.function(p.value)) p.value <- p.value(observed, residuals, expected, df)
    v <- if(p.value < (1-level)) my.v[1] else my.v[2]
  }

  ## set up function for interpolation of saturation
  if(!is.function(interpolate)) {
    col.bins <- sort(interpolate)
    interpolate <- stepfun(col.bins,  seq(my.s[2], my.s[1], length = length(col.bins) + 1))
    col.bins <- sort(unique(c(col.bins, 0, -col.bins)))
  } else {
    col.bins <- NULL
  }

  ## store color and lty information for legend
  legend <- NULL
  if(!is.null(col.bins)) {
    res2 <- col.bins
    res2 <- c(head(res2, 1) - 1, res2[-1] - diff(res2)/2, tail(res2, 1) + 1)
    legend.col <- hsv(ifelse(res2 > 0, my.h[1], my.h[2]),
                      pmax(pmin(interpolate(abs(res2)), 1), 0),
		      v, ...)
    lty.bins <- 0
    legend.lty <- lty[2:1]
    legend <- list(col = legend.col, col.bins = col.bins,
                   lty = legend.lty, lty.bins = lty.bins)
  }

  ## set up function that computes color/lty from residuals
  rval <- function(x) {
    res <- as.vector(x)

    fill <- hsv(ifelse(res > 0, my.h[1], my.h[2]),
                pmax(pmin(interpolate(abs(res)), 1), 0),
	        v, ...)
    dim(fill) <- dim(x)

    col <- rep(line_col, length.out = length(res))
    if(!is.null(eps)) {
      eps <- abs(eps)
      col[res > eps] <- hsv(my.h[1], 1, v, ...)
      col[res < -eps] <- hsv(my.h[2], 1, v, ...)
    }
    dim(col) <- dim(x)

	# line type should be solid if abs(resid) < eps
	ltytmp <- ifelse(x > 0, lty[1], lty[2])
	if(!is.null(eps))
		ltytmp[abs(x) < abs(eps)] <- lty[1]
	dim(ltytmp) <- dim(x)

	return(structure(list(col = col, fill = fill, lty = ltytmp), class = "gpar"))
}
  attr(rval, "legend") <- legend
  attr(rval, "p.value") <- p.value
  return(rval)
}
class(shading_hsv) <- "grapcon_generator"


shading_hcl <- function(observed, residuals = NULL, expected = NULL, df = NULL,
  h = NULL, c = NULL, l = NULL,
  interpolate = c(2, 4), lty = 1, eps = NULL, line_col = "black",
  p.value = NULL, level = 0.95, ...)
{
  ## set defaults
  if(is.null(h)) h <- c(260, 0)
  if(is.null(c)) c <- c(100, 20)
  if(is.null(l)) l <- c(90, 50)

  ## get h/c/l and lty
  my.h <- rep(h, length.out = 2)  ## positive and negative hue
  my.c <- rep(c, length.out = 2)  ## significant and non-significant maximum chroma
  my.l <- rep(l, length.out = 2)  ## maximum and minimum luminance
  lty <- rep(lty, length.out = 2) ## positive and negative lty

  ## model fitting (if necessary)
  if(is.null(expected) && !is.null(residuals)) stop("residuals without expected values specified")
  if(!is.null(expected) && is.null(df) && is.null(p.value)) {
    warning("no default inference available without degrees of freedom")
    p.value <- NA
  }
  if(is.null(expected) && !is.null(observed)) {
    expected <- loglin(observed, 1:length(dim(observed)), fit = TRUE, print = FALSE)
    df <- expected$df
    expected <- expected$fit
  }
  if(is.null(residuals) && !is.null(observed)) residuals <- (observed - expected)/sqrt(expected)

  ## conduct significance test (if specified)
  if(is.null(p.value)) p.value <- function(observed, residuals, expected, df)
    pchisq(sum(as.vector(residuals)^2), df, lower.tail = FALSE)
  if(!is.function(p.value) && is.na(p.value)) {
    max.c <- my.c[1]
    p.value <- NULL
  } else {
    if(is.function(p.value)) p.value <- p.value(observed, residuals, expected, df)
    max.c <- ifelse(p.value < (1-level), my.c[1], my.c[2])
  }

  ## set up function for interpolation of saturation
  if(!is.function(interpolate)) {
    col.bins <- sort(interpolate)
    interpolate <- stepfun(col.bins,  seq(0, 1, length = length(col.bins) + 1))
    col.bins <- sort(unique(c(col.bins, 0, -col.bins)))
  } else {
    col.bins <- NULL
  }

  ## store color and lty information for legend
  legend <- NULL
  if(!is.null(col.bins)) {
    res2 <- col.bins
    res2 <- c(head(res2, 1) - 1, res2[-1] - diff(res2)/2, tail(res2, 1) + 1)
    legend.col <- hcl2hex(ifelse(res2 > 0, my.h[1], my.h[2]),
                      max.c * pmax(pmin(interpolate(abs(res2)), 1), 0),
	              my.l[1] + diff(my.l) * pmax(pmin(interpolate(abs(res2)), 1), 0),
		      ...)
    lty.bins <- 0
    legend.lty <- lty[2:1]
    legend <- list(col = legend.col, col.bins = col.bins,
                   lty = legend.lty, lty.bins = lty.bins)
  }

  ## set up function that computes color/lty from residuals
  rval <- function(x) {
    res <- as.vector(x)

    fill <- hcl2hex(ifelse(res > 0, my.h[1], my.h[2]),
                max.c * pmax(pmin(interpolate(abs(res)), 1), 0),
	        my.l[1] + diff(my.l) * pmax(pmin(interpolate(abs(res)), 1), 0),
	        ...)
    dim(fill) <- dim(x)

    col <- rep(line_col, length.out = length(res))
    if(!is.null(eps)) {
      eps <- abs(eps)
      col[res > eps] <- hcl2hex(my.h[1], max.c, my.l[2], ...)
      col[res < -eps] <- hcl2hex(my.h[2], max.c, my.l[2], ...)
    }
    dim(col) <- dim(x)

    ltytmp <- ifelse(x > 0, lty[1], lty[2])
    if(!is.null(eps))
        ltytmp[abs(x) < abs(eps)] <- lty[1]
    dim(ltytmp) <- dim(x)

    return(structure(list(col = col, fill = fill, lty = ltytmp), class = "gpar"))
  }
  attr(rval, "legend") <- legend
  attr(rval, "p.value") <- p.value
  return(rval)
}
class(shading_hcl) <- "grapcon_generator"

shading_Friendly <- function(observed = NULL, residuals = NULL, expected = NULL, df = NULL,
  h = c(2/3, 0), lty = 1:2, interpolate = c(2, 4), eps = 0.01, line_col = "black", ...)
{
  shading_hsv(observed = NULL, residuals = NULL, expected = NULL, df = NULL,
              h = h, v = 1, lty = lty, interpolate = interpolate,
	      eps = eps, line_col = line_col, p.value = NA, ...)
}
class(shading_Friendly) <- "grapcon_generator"

shading_Friendly2 <- function(observed = NULL, residuals = NULL, expected = NULL, df = NULL, lty = 1:2, interpolate = c(2, 4), eps = 0.01, line_col = "black", ...)
{
  shading_hcl(observed = NULL, residuals = NULL, expected = NULL, df = NULL,
              lty = lty, interpolate = interpolate,
	      eps = eps, line_col = line_col, p.value = NA, ...)
}
class(shading_Friendly2) <- "grapcon_generator"

shading_sieve <-
    function(observed = NULL, residuals = NULL, expected = NULL, df = NULL,
             h = c(260, 0), lty = 1:2, interpolate = c(2, 4), eps = 0.01,
             line_col = "black", ...)
{
    shading_hcl(observed = NULL, residuals = NULL, expected = NULL,
                df = NULL, h = h, c = 100, l = 50, lty = lty,
                interpolate = interpolate,
                eps = eps, line_col = line_col, p.value = NA, ...)
}
class(shading_sieve) <- "grapcon_generator"

shading_max <- function(observed = NULL, residuals = NULL, expected = NULL, df = NULL,
  h = NULL, c = NULL, l = NULL, lty = 1, eps = NULL, line_col = "black", level = c(0.9, 0.99), n = 1000, ...)
{
  stopifnot(length(dim(observed)) == 2)

  ## set defaults
  if(is.null(h)) h <- c(260, 0)
  if(is.null(c)) c <- c(100, 20)
  if(is.null(l)) l <- c(90, 50)

  obs.test <- coindep_test(observed, n = n)
  col.bins <- obs.test$qdist(sort(level))
  rval <- shading_hcl(observed = NULL, residuals = NULL, expected = NULL, df = NULL,
                        h = h, c = c, l = l, interpolate = col.bins, lty = lty,
			eps = eps, line_col = line_col, p.value = obs.test$p.value, ...)
  return(rval)
}
class(shading_max) <- "grapcon_generator"

shading_binary <- function(observed = NULL, residuals = NULL, expected = NULL, df = NULL,
  col = NULL)
{
  ## check col argument
  if(is.null(col)) col <- hcl2hex(c(260, 0), 50, 70)
  col <- rep(col, length.out = 2)

  ## store color information for legend
  legend <- list(col = col[2:1], col.bins = 0, lty = NULL, lty.bins = NULL)

  ## set up function that computes color/lty from residuals
  rval <- function(x)
    gpar(fill = ifelse(x > 0, col[1], col[2]))

  ## add meta information for legend
  attr(rval, "legend") <- legend
  attr(rval, "p.value") <- NULL

  rval
}
class(shading_binary) <- "grapcon_generator"

shading_Marimekko <-
function(x, fill = NULL, byrow = FALSE)
{
    if (is.null(fill)) fill <- colorspace::rainbow_hcl
    d <- dim(x)
    l1 <- if (length(d) > 1L) d[2] else d
    l2 <- if (length(d) > 1L) d[1] else 1
    if (is.function(fill)) fill <- fill(l1)
    fill <- if (byrow) rep(fill, l2) else rep(fill, each = l2)
    gpar(col = NA, lty = "solid",
         fill = array(fill, dim = d))
}

shading_diagonal <-
function(x, fill = NULL)
{
    if (is.null(fill)) fill <- colorspace::rainbow_hcl
    d <- dim(x)
    if (length(d) < 1L)
        stop("Need matrix or array!")
    if (d[1] != d[2])
        stop("First two dimensions need to be of same length!")
    if (is.function(fill)) fill <- fill(d[1])
    tp = toeplitz(seq_len(d[1]))
    gpar(col = NA, lty = "solid",
         fill = array(rep(fill[tp],  d[1]), dim = d))
}
