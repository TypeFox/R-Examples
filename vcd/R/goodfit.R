goodfit <- function(x, type = c("poisson", "binomial", "nbinomial"),
                    method = c("ML", "MinChisq"), par = NULL)
{
    if(is.vector(x)) {
        x <- table(x)
    }
    if(is.table(x)) {
        if(length(dim(x)) > 1) stop ("x must be a 1-way table")
        freq <- as.vector(x)
        count <- as.numeric(names(x))
    } else {
          if(!(!is.null(ncol(x)) && ncol(x) == 2))
              stop("x must be a 2-column matrix or data.frame")
          freq <- as.vector(x[,1])
          count <- as.vector(x[,2])
      }

    ## fill-in possibly missing cells
    nfreq <- rep(0, max(count) + 1)
    nfreq[count + 1] <- freq
    freq <- nfreq
    count <- 0:max(count)
    n <- length(count)

    ## starting value for degrees of freedom
    df <- -1

    type <- match.arg(type)
    method <- match.arg(method)

    switch(type,

           "poisson" = {

               if(!is.null(par)) {
                   if(!is.list(par))
                       stop("`par' must be a named list")
                   if(names(par) != "lambda")
                       stop("`par' must specify `lambda'")

                   par <- par$lambda
                   method <- "fixed"
               }
               else if(method == "ML") {
                   df <- df - 1
                   par <- weighted.mean(count,freq)
               }
               else if(method == "MinChisq") {
                   df <- df - 1

                   chi2 <- function(x)
                   {
                   p.hat <- diff(c(0, ppois(count[-n], lambda = x), 1))
                   expected <- sum(freq) * p.hat
                   sum((freq - expected)^2/expected)
                   }

                   par <- optimize(chi2, range(count))$minimum
               }
               par <- list(lambda = par)
               p.hat <- dpois(count, lambda = par$lambda)
           },

           "binomial" = {
               size <- par$size
               if(is.null(size)) {
                   size <- max(count)
                   warning("size was not given, taken as maximum count")
               }

               if(size > max(count)) {
                   nfreq <- rep(0, size + 1)
                   nfreq[count + 1] <- freq
                   freq <- nfreq
                   count <- 0:size
                   n <- length(count)
               }

               if(!is.null(par$prob)) {
                   if(!is.list(par))
                       stop("`par' must be a named list and specify `prob'")
                   par <- par$prob
                   method <- "fixed"
               }
               else if(method == "ML") {
                   df <- df - 1
                   par <- weighted.mean(count/size, freq)
               }
               else if(method == "MinChisq") {
                   df <- df - 1

                   chi2 <- function(x)
                   {
                   p.hat <- diff(c(0, pbinom(count[-n], prob = x, size = size), 1))
                   expected <- sum(freq) * p.hat
                   sum((freq - expected)^2/expected)
                   }

                   par <- optimize(chi2, c(0,1))$minimum
               }
               par <- list(prob = par, size = size)
               p.hat <- dbinom(count, prob = par$prob, size = par$size)
           },


           "nbinomial" = {

               if(!is.null(par)) {
                   if(!is.list(par)) stop("`par' must be a named list")
                   if(!(isTRUE(all.equal(names(par), "size")) | isTRUE(all.equal(sort(names(par)), c("prob", "size")))))
                       stop("`par' must specify `size' and possibly `prob'")
                   if(!is.null(par$prob)) method <- "fixed"
               }

               switch(method,

                      "ML" = {
                          if(is.null(par$size)) {
                              df <- df - 2
                              par <- fitdistr(rep(count, freq), "negative binomial")$estimate
                              par <- par[1]/c(1, sum(par))
                          } else {
                                df <- df - 1
                                method <- c("ML", "with size fixed")

                                size <- par$size
                                xbar <- weighted.mean(count,freq)
                                par <- c(size, size/(xbar+size))
                            }
                      },

                      "MinChisq" = {
                          if(is.null(par$size)) {
                              df <- df - 2

                              ## MM
                              xbar <- weighted.mean(count,freq)
                              s2 <- var(rep(count,freq))
                              p <- xbar / s2
                              size <- xbar^2/(s2 - xbar)
                              par1 <- c(size, p)

                              ## minChisq
                              chi2 <- function(x)
                              {
                              p.hat <- diff(c(0, pnbinom(count[-n], size = x[1], prob = x[2]), 1))
                              expected <- sum(freq) * p.hat
                              sum((freq - expected)^2/expected)
                              }

                              par <- optim(par1, chi2)$par
                          } else {
                                df <- df - 1
                                method <- c("MinChisq", "with size fixed")

                                chi2 <- function(x)
                                {
                                p.hat <- diff(c(0, pnbinom(count[-n], size = par$size, prob = x), 1))
                                expected <- sum(freq) * p.hat
                                sum((freq - expected)^2/expected)
}
                                par <- c(par$size, optimize(chi2, c(0, 1))$minimum)
                            }
                      },

                      "fixed" = {
                          par <- c(par$size, par$prob)
                      })

               par <- list(size = par[1], prob = par[2])
               p.hat <- dnbinom(count, size = par$size, prob = par$prob)
           })

    expected <- sum(freq) * p.hat

    df <- switch(method[1],
                 "MinChisq" = { length(freq) + df },
                 "ML" = { sum(freq > 0) + df },
                 "fixed" = { c(length(freq), sum(freq > 0)) + df }
                 )

    structure(list(observed = freq,
                   count = count,
                   fitted = expected,
                   type = type,
                   method = method,
                   df = df,
                   par = par),
              class = "goodfit")
}

# does this need a residuals_type arg?
print.goodfit <- function(x, residuals_type = c("pearson", "deviance", "raw"), ...)
{
    residuals_type <- match.arg(residuals_type)
    cat(paste("\nObserved and fitted values for", x$type, "distribution\n"))
    if(x$method[1] == "fixed")
      cat("with fixed parameters \n\n")
    else
      cat(paste("with parameters estimated by `", paste(x$method, collapse = " "), "' \n\n", sep = ""))

    resids <- residuals(x, type = residuals_type)
    RVAL <- cbind(x$count, x$observed, x$fitted, resids)
    colnames(RVAL) <- c("count", "observed", "fitted",
                        paste(residuals_type, "residual"))
    rownames(RVAL) <- rep("", nrow(RVAL))
    print(RVAL, ...)
    invisible(x)
}

summary.goodfit <- function(object, ...)
{
    df <- object$df
    obsrvd <- object$observed
    count <- object$count
    expctd <- fitted(object)

    G2 <- sum(ifelse(obsrvd == 0, 0, obsrvd * log(obsrvd/expctd))) * 2

    n <- length(obsrvd)
    pfun <- switch(object$type,
                   poisson = "ppois",
                   binomial = "pbinom",
                   nbinomial = "pnbinom")
    p.hat <- diff(c(0, do.call(pfun, c(list(q = count[-n]), object$par)), 1))
    expctd <- p.hat * sum(obsrvd)
    X2 <- sum((obsrvd - expctd)^2 / expctd)

    names(G2) <- "Likelihood Ratio"
    names(X2) <- "Pearson"
    if(any(expctd < 5) & object$method[1] != "ML")
        warning("Chi-squared approximation may be incorrect")

    RVAL <- switch(object$method[1],
                   ML = G2,
                   MinChisq = X2,
                   fixed = c(X2, G2)
                   )

    RVAL <- cbind(RVAL, df, pchisq(RVAL, df = df, lower.tail = FALSE))
    colnames(RVAL) <- c("X^2", "df", "P(> X^2)")

    cat(paste("\n\t Goodness-of-fit test for", object$type, "distribution\n\n"))
    print(RVAL, ...)
    invisible(RVAL)
}

plot.goodfit <- function(x, ...)
{
    rootogram(x, ...)
}

fitted.goodfit <- function(object, ...)
{
    object$fitted
}

residuals.goodfit <- function(object, type = c("pearson", "deviance", "raw"),  ...)
{
    obsrvd <- object$observed
    expctd <- fitted(object)
    count <- object$count
    n <- length(obsrvd)
    pfun <- switch(object$type,
                   poisson = "ppois",
                   binomial = "pbinom",
                   nbinomial = "pnbinom")
    p.hat <- diff(c(0, do.call(pfun, c(list(q = count[-n]), object$par)), 1))
    expctd <- p.hat * sum(obsrvd)

    res <- switch(match.arg(type),
                  pearson = (obsrvd - expctd) / sqrt(expctd),
                  deviance = ifelse(obsrvd == 0, 0,
                                    obsrvd * log(obsrvd / expctd)),
                  obsrvd - expctd)

    return(res)
}

predict.goodfit <- function(object, newcount = NULL, type = c("response", "prob"), ...)
{
    if(is.null(newcount)) newcount <- object$count
    type <- match.arg(type)

    densfun <- switch(object$type,
                      poisson = "dpois",
                      binomial = "dbinom",
                      nbinomial = "dnbinom")

    RVAL <- do.call(densfun, c(list(x = newcount), object$par))
    if (type == "response") RVAL <- RVAL * sum(object$observed)
    return(RVAL)
}

