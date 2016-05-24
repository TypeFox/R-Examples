# Calculate the 5-point stencil
# point f'(x.i)
first.midpoint <- function(y, h)
  1/12/h * (y[1] - 8 * y[2] + 8 * y[4] - y[5])

# point f'(x.0)
first.beginpoint0 <- function(y, h)
  1/12/h * (-25 * y[1] + 48 * y[2] - 36 * y[3] + 16 * y[4] - 3 * y[5])

# point f'(x.1)
first.beginpoint1 <- function(y, h)
  1/12/h * (-3 * y[1] - 10 * y[2] + 18 * y[3] - 6 * y[4] + y[5])

# point f'(x.{n-1})
first.endpoint1 <- function(y, h)
  - first.beginpoint1(rev(y), h)

# point f'(x.n)
first.endpoint0 <- function(y, h)
  - first.beginpoint0(rev(y), h)

# point f''(x.i)
sec.midpoint <- function(y, h)
  1/12/h^2 * (- y[1] + 16 * y[2] - 30 * y[3]  +  16 * y[4] - y[5])

# point f''(x.0)
sec.beginpoint0 <- function(y, h)
  1/12/h^2 * (35 * y[1] - 104 * y[2] + 114 * y[3]  -  56 * y[4] + 11 * y[5])

# point f''(x.1)
sec.beginpoint1 <- function(y, h)
  1/12/h^2 * (11 * y[1] - 20 * y[2] + 6 * y[3]  +  4 * y[4] - y[5])

# point f''(x.{n-1})
sec.endpoint1 <- function(y, h)
  sec.beginpoint1(rev(y), h)

# point f''(x.n)
sec.endpoint0 <- function(y, h)
  sec.beginpoint0(rev(y), h)

# Defintion of the inder function
inder <- function(x, y, Nip = 4, logy = FALSE, smooth.method = "spline") {
  # Test validity ot the input data
  testxy(x, y)
  
  # Test meaningfulness of the spline interpolation and give a warning in case
  # any violation
  if (Nip < 1) 
    stop("Use Nip equal or larger to 1")
  
  if (Nip > 10) 
    warning("Nip larger than 10 may case over-fitting")
  
  # Set the smooth method for inder
  
  if(is.null(smooth.method)) {
    tmp.xy <- data.frame(x = x, y = y)
    smooth.method <- "no smoothing"
  } else {
    if (smooth.method == "spline") {
      tmp.xy <- spline(x, y, n = Nip * length(x))
    } 
    
    if (smooth.method == "supsmu") {
      tmp.xy <- supsmu(x, y, span = "cv")
    } 
  }
  
  x <- tmp.xy[["x"]]
  y <- tmp.xy[["y"]]
  
  if (logy == TRUE) y <- log10(tmp.xy[["y"]])
  
  # calculate h, naive approach
  h <- vapply(2L:length(x), function(i) x[i] - x[i - 1], 0)
  # instead of zero, in statement should be the minina
  if (var(h) > .Machine[["double.eps"]]) 
    warning("Points are not equidistant. The results of interpolation 
	      could be not correct.")
  
  h <- h[1]
  # calculate midpoints
  first.der <- c(first.beginpoint0(y[1:5], h),
                 first.beginpoint1(y[1:5], h),
                 vapply(3L:(length(y) - 2), function(i) first.midpoint(y[(i - 2):(i + 2)], h), 0),
                 first.endpoint1(y[(length(y) - 5):length(y)], h),
                 first.endpoint0(y[(length(y) - 5):length(y)], h))
  
  sec.der <- c(sec.beginpoint0(y[1:5], h),
               sec.beginpoint1(y[1:5], h),
               vapply(3L:(length(y) - 2), function(i) sec.midpoint(y[(i - 2):(i + 2)], h), 0),
               sec.endpoint1(y[(length(y) - 5):length(y)], h),
               sec.endpoint0(y[(length(y) - 5):length(y)], h))
  
  dat <- cbind(x, y, first.der, sec.der)
  colnames(dat) <- c("x", "y", "d1y", "d2y")
  new("der", '.Data' = dat, 'method' = smooth.method)
}

setGeneric("inder")

setMethod("inder", signature(x = "data.frame", y="missing"), 
          function(x, y, Nip = 4, logy = FALSE, smooth.method = "spline") { 
            if (ncol(x) != 2) 
              stop("'x' must have two columns.")
            inder(x[, 1], x[, 2], Nip = Nip, logy = logy, 
                  smooth.method = smooth.method)
          })

setMethod("inder", signature(x = "matrix", y = "missing"), 
          function(x, y, Nip = 4, logy = FALSE, smooth.method = "spline") { 
            if (ncol(x) != 2) 
              stop("'x' must have two columns.")
            inder(x[, 1], x[, 2], Nip = Nip, logy = logy, 
                  smooth.method = smooth.method)
          })
