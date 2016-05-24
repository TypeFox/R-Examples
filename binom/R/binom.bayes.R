binom.bayes <- function(x, n,
                        conf.level = 0.95,
                        type = c("highest", "central"),
                        prior.shape1 = 0.5,
                        prior.shape2 = 0.5, # Jeffrey's prior
                        tol = .Machine$double.eps^.5,
                        maxit = 1000, ...) {
  if(prior.shape1 <= 0 || prior.shape2 <= 0)
    stop("priors must be strictly positive.")
  if((length(x) != length(n)) && length(x) == 1)
    x <- rep(x, length(n))
  if((length(x) != length(n)) && length(n) == 1)
    n <- rep(n, length(x))
  ends <- x == 0 | x == n
  alpha <- rep(1 - conf.level, length(x))
  alpha[!ends] <- alpha[!ends] * 0.5
  a <- x + prior.shape1
  b <- n - x + prior.shape2
  p <- a/(a + b)
  lcl <- qbeta(alpha, a, b)
  ucl <- qbeta(1 - alpha, a, b)
  type <- match.arg(type)
  error <- vector("logical", length(lcl))
  if(any(!ends) && (type == "highest")) {
    ci <- .C("binom_bayes",
             as.integer(x[!ends]),
             as.integer(n[!ends]),
             as.double(a[!ends]),
             as.double(b[!ends]),
             as.double(alpha[!ends]),
             lcl = as.double(lcl[!ends]),
             ucl = as.double(ucl[!ends]),
             as.integer(sum(!ends)),
             maxit = as.integer(maxit),
             tol = as.double(tol),
             error = error[!ends],
             PACKAGE = "binom")
    lcl[!ends] <- ci$lcl
    ucl[!ends] <- ci$ucl
    error[!ends] <- as.logical(ci$error)
    if(any(ci$error)) {
      nerr <- sum(ci$error)
      msg1 <- sprintf("%s confidence interval%s ", nerr, if(nerr > 1) "s" else "")
      msg2 <- "failed to converge (marked by '*').\n"
      msg3 <- "  Try changing 'tol' to a different value."
      warning(msg1, msg2, msg3)
    }
  }
  lcl[x == 0] <- 0
  ucl[x == n] <- 1
  sig <- pbeta(lcl, a, b) + 1 - pbeta(ucl, a, b)
  res <- data.frame(method = "bayes",
                    x = x,
                    n = n,
                    shape1 = a,
                    shape2 = b,
                    mean = p,
                    lower = lcl,
                    upper = ucl,
                    sig = sig)
  if(any(error))
    res$method <- factor(sprintf("bayes%s", ifelse(ci$error, "*", "")))
  attr(res, "conf.level") <- conf.level
  attr(res, "type") <- type
  res
}

rbind.zero <- function(data, which.x = c("max", "min", "both"), row.only = FALSE) {
  if (nrow(data) == 0) return(data)
  which.x <- match.arg(which.x)
  if (which.x == "both") {
    data.max <- rbind.zero(data, which.x = "max", row.only = TRUE)
    data.min <- rbind.zero(data, which.x = "min", row.only = TRUE)
    data.x <- rbind(data.max, data.min)
  } else {
    which.fn <- if (which.x == "max") which.max else which.min
    data.x <- data[which.fn(data$xx), ]
    data.x$yy <- 0
  }
  if (row.only) data.x else rbind(data, data.x)
}

build.density.data <- function(i, x, bayes) {
  y <- dbeta(x, bayes$shape1[i], bayes$shape2[i])
  b <- bayes[rep(i, length(x)), ]
  data <- data.frame(xx = x, yy = y, b)
  data <- data[is.finite(data$yy), ]
  all.labels <- sprintf("x = %d, n = %d, conf.level = %0.2f",
                        bayes$x, bayes$n, 1 - bayes$sig)
  all.labels <- reorder(all.labels, bayes$mean)
  data$label <- all.labels[i]
  row.names(data) <- NULL
  data.lower <- data[data$xx < data$lower, ]
  data.upper <- data[data$xx > data$upper, ]
  data.central <- data[data$xx >= data$lower & data$xx <= data$upper, ]
  data.lower <- rbind.zero(data.lower, which.x = "max")
  data.upper <- rbind.zero(data.upper, which.x = "min")
  data.central <- rbind.zero(data.central, which.x = "both")
  list(lower = data.lower, upper = data.upper, central = data.central)
}

binom.bayes.densityplot <- function(bayes,
                                    npoints = 500,
                                    fill.central = "lightgray",
                                    fill.lower = "steelblue",
                                    fill.upper = fill.lower,
                                    alpha = 0.8, ...) {
  stopifnot(require(ggplot2))
  x <- seq(0, 1, length.out = npoints)
  datalist <- lapply(1:nrow(bayes), build.density.data, x = x, bayes = bayes)
  data.lower <- do.call(rbind, lapply(datalist, "[[", "lower"))
  data.upper <- do.call(rbind, lapply(datalist, "[[", "upper"))
  data.central <- do.call(rbind, lapply(datalist, "[[", "central"))
  gg <- ggplot(data.central, aes_string(x = "xx", y = "yy"))
  gg <- gg + geom_polygon(fill = fill.central, alpha = alpha)
  if (nrow(data.lower) > 0) {
    gg <- gg + geom_polygon(data = data.lower, fill = fill.lower, alpha = alpha)
  }
  if (nrow(data.upper) > 0) {
    gg <- gg + geom_polygon(data = data.upper, fill = fill.upper, alpha = alpha)
  }
  gg <- gg + facet_wrap(~ label)
  gg <- gg + xlim(c(0, 1)) + xlab(NULL)
  gg <- gg + ylim(c(0, max(data.central$y))) + ylab("Beta Density")
  gg <- gg + theme_bw()
  gg
}
