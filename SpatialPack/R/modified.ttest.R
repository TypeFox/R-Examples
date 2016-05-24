modified.ttest <-
function(x, y, coords, nclass = 13)
{
  ## validating arguments
  if (length(x) != length(y))
    stop("'x' and 'y' must have the same length")
  if (!is.numeric(x)) stop("'x' must be a numeric vector")
  if (!is.numeric(y)) stop("'y' must be a numeric vector")
  ## in order to remove all NAs
  OK <- complete.cases(x, y)
  x <- x[OK]
  y <- y[OK]
  n <- length(x)
  corr <- cor(cbind(x, y))
  dnames <- colnames(corr)
  corr <- corr[2,1]

  ## extract coordinates, is assumed that the variables are in the appropiate order
  coords <- as.matrix(coords)
  p <- ncol(coords)
  if (p < 2) stop("'coords' must be a matrix with two columns")
  if (p > 2) warning("only the first two columns of 'coords' are considered")
  p <- 2 # only implemented for this case!
  xpos <- coords[,1]
  ypos <- coords[,2]
  cnames <- colnames(coords)[1:p]

  ## some definitions
  ndist  <- n * (n - 1) / 2
  if (is.null(nclass))
    nclass = as.integer(1.5 + 3.3 * log10(ndist))
  dims <- c(n, p, nclass)

  ## call routine
  now <- proc.time()
  z <- .C("modified_ttest",
          x = as.double(x),
          y = as.double(y),
          xpos = as.double(xpos),
          ypos = as.double(ypos),
          dims = as.integer(dims),
          corr = as.double(corr),
          upper.bounds = double(nclass),
          card = double(nclass),
          imoran = double(nclass * p),
          stats = double(3))
  speed <- proc.time() - now

  ## creating output object
  o <- list(corr = z$corr, Fstat = z$stats[1], dof = z$stats[2], p.value = z$stats[3])
  o$dims <- dims
  o$upper.bounds <- z$upper.bounds
  o$card <- z$card
  o$imoran <- matrix(z$imoran, ncol = p)
  colnames(o$imoran) <- dnames
  o$imoran <- as.data.frame(o$imoran)
  o$data.names <- dnames
  o$coords.names <- cnames
  o$speed <- speed
  class(o) <- "mod.ttest"
  return(o)
}

print.mod.ttest <- function(x, digits = 4, ...)
{
  cat("\n")
  cat("Corrected Pearson's correlation for spatial autocorrelation\n")
  cat("\n")
  dnames <- paste(x$data.names, collapse = " and ", sep = "")
  coords <- paste(x$coords.names, collapse = " and ", sep = "")
  cat("data:", paste(c(dnames, ";"), sep = ""), "coordinates:", coords, "\n")
  cat("F-statistic:", format(round(x$Fstat, digits = digits)), "on 1 and",
      format(round(x$dof, digits = digits)), "DF, p-value:",
      format(round(x$p.value, digits = digits)), "\n")
  cat("alternative hypothesis: true autocorrelation is not equal to 0\n")
  cat("sample correlation:", format(round(x$corr, digits = digits)))
  cat("\n")
  invisible(x)
}

summary.mod.ttest <- function(object, ...)
{
  z <- object
  coef <- cbind(z$upper.bounds, z$card, as.matrix(z$imoran))
  nclass <- z$dims[3]
  dimnames(coef) <- list(1:nclass, c("Upper Bounds", "Cardinality", "Moran:x", "Moran:y"))
  ans <- z[c("corr", "Fstat", "dof", "p.value")]
  ans$data.names <- z$data.names
  ans$coords.names <- z$coords.names
  ans$coef <- coef
  class(ans) <- "summary.mod.ttest"
  ans
}

print.summary.mod.ttest <- function(x, digits = 4, ...)
{
  print.mod.ttest(x)
  cat("\n")
  print(x$coef, digits = digits)
  invisible(x)
}
