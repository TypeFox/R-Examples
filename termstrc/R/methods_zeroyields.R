### Nelson/Siegel and Svensson spot curve estimaton from zero yields 

zeroyields <- function(maturities, yields, dates)
  {
    zy <- list(maturities = maturities, yields = yields, dates = dates)
    class(zy) <- "zeroyields"
    zy
  }

print.zeroyields <- function(x, ...)
  {
    cat("This is a data set of zero-coupon yields.\n")
    cat(paste("Maturities range from", min(x$maturities), "to", max(x$maturities),"years.\n"))
    cat(paste("There are",nrow(x$yields), "observations between",x$dates[1], "and",x$dates[length(x$dates)],".\n"))
  }

summary.zeroyields <- function(object, ...)
  {
    print(summary(object$yields))
  }

plot.zeroyields <- function(x,...)
  {
    Z <- as.matrix(x$yields)
    X <- 1:nrow(x$yields)
    Y <- x$maturities

    open3d()
    persp3d(X, Y, Z, col = "green3", box = FALSE,xlab = "Time", ylab = "Maturities (years)", zlab = "Zero-yields (%)")
  }

