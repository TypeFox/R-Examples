`plotLtt` <-
function(x)
{
  if (!is.numeric(x)) stop("object x not of class 'numeric'")
  x <- rev(sort(x))
  x <- c(x, 0)
  n <- length(x)
  z <- max(x)-(x)

  vn <- 1:n
  plot(log(vn+1)~z, xlab = "Time From Basal Divergence", ylab = "Log Lineages",
    main = "Log-Lineages Through Time")
  print(log(vn+1))
}

