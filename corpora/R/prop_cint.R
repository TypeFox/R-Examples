prop.cint <- function(k, n, method=c("binomial", "z.score"), correct=TRUE,
                      conf.level=0.95, alternative=c("two.sided", "less", "greater")) {
  method <- match.arg(method)
  alternative <- match.arg(alternative)
  if (any(k < 0) || any(k > n) || any(n < 1)) stop("arguments must be integer vectors with 0 <= k <= n")
  if (any(conf.level <= 0) || any(conf.level > 1)) stop("conf.level must be in range [0,1]")

  l <- max(length(k), length(n), length(conf.level)) # ensure that all vectors have the same length
  if (length(k) < l) k <- rep(k, length.out=l)
  if (length(n) < l) n <- rep(n, length.out=l)
  if (length(conf.level) < l) conf.level <- rep(conf.level, length.out=l)

  if (method == "binomial") {
    ## compute binomial confidence interval (using incomplete Beta function)
    alpha <- if (alternative == "two.sided") (1 - conf.level) / 2 else (1 - conf.level)
    lower <- qbeta(alpha, k, n - k + 1)
    upper <- qbeta(1 - alpha, k + 1, n - k)
    cint <- switch(alternative,
                   two.sided = data.frame(lower = lower, upper = upper),
                   less      = data.frame(lower = 0,     upper = upper),
                   greater   = data.frame(lower = lower, upper = 1))
  } else {
    ## compute z-score confidence interval (by solving quadratic z-test equation for p)
    alpha <- if (alternative == "two.sided") (1 - conf.level) / 2 else (1 - conf.level)
    z <- qnorm(alpha, lower.tail=FALSE) # z-score corresponding to desired confidence level
    yates <- if (correct) 0.5 else 0.0  # whether to apply Yates' correction
    
    k.star <- k - yates                 # lower boundary of confidence interval (solve implicit equation for z-score test)
    k.star <- pmax(0, k.star)           # Yates' correction cannot be satisfied at boundary of valid range for k
    A <- n + z^2                        # coefficients of quadratic equation that has to be solved
    B <- -2 * k.star - z^2
    C <- k.star^2 / n
    lower <- solve.quadratic(A, B, C, nan.lower=0)$lower

    k.star <- k + yates                 # upper boundary of confidence interval
    k.star <- pmin(n, k.star)
    A <- n + z^2
    B <- -2 * k.star - z^2
    C <- k.star^2 / n
    upper <- solve.quadratic(A, B, C, nan.upper=1)$upper

    cint <- switch(alternative,
                   two.sided = data.frame(lower = lower,    upper = upper),
                   less      = data.frame(lower = rep(0,l), upper = upper),
                   greater   = data.frame(lower = lower,    upper = rep(1,l)))
  }

  cint
}
