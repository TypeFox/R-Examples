sbdiv <-
  function(X, f, theta = c("Shannon", "Simpson"), type = c("Dunnett", "Tukey", "Sequen", "AVE", "Changepoint", "Williams", "Marcus", "McDermott", "UmbrellaWilliams", "GrandMean"), cmat = NULL, method = c("WYht", "tsht", "rpht", "asht"), conf.level = 0.95, alternative = c("two.sided", "less", "greater"), R = 2000, base = 1,...)
  {
    args <- list(...)
    theta <- match.arg(theta)
    type <- match.arg(type)
    method <- match.arg(method)
    alternative <- match.arg(alternative)
    n <- unlist(lapply(split(x = X, f = f), length))
    if (length(conf.level) != 1) {
      stop("argument conf.level should be a single numeric value")
    }
    if (conf.level <= 0.5 | conf.level >= 1) {
      stop("argument conf.level should be a single numeric value between 0.5 and 1")
    }
    if (!is.numeric(conf.level)) {
      stop("argument conf.level should be a single numeric value")
    }
    if (!is.data.frame(X)) {
      stop("X must be of an object of class 'data.frame'")
    }
    if (!is.factor(f)) {
      f <- as.factor(f)
    }
    k <- length(levels(f))
    if (k <= 1) {
      stop("The factor variable f should have at least 2 levels to be compared")
    }
    if (is.null(cmat))
      {
        cmat <- contrMat(n = n, type = type, base = base)
      }
    else {
      if (ncol(cmat) != k) {
        stop("Number of columns in cmat should be the same as the number of levels in f")
      }
    }
    switch(theta,
           "Shannon" = {
             switch(method,
                    "WYht" = {
                      output <- WYht(X = X, f = f, theta = estShannonWY, cmat = cmat, conf.level = conf.level, alternative = alternative, R = R, args = args)
                    },
                    "tsht" = {
                      output <- tsht(X = X, f = f, theta = estShannonf, cmat = cmat, conf.level = conf.level, alternative = alternative, R = R, args = args)
                    },
                    "rpht" = {
                      output <- rpht(X = X, f = f, theta = estShannonf, cmat = cmat, conf.level = conf.level, alternative = alternative, R = R, args = args)
                    },
                    "asht" = {
                      output <- asht(X = X, f = f, theta = estShannonf, cmat = cmat, conf.level = conf.level, alternative = alternative, args = args)
                    }
                    )
           },
           "Simpson" = {
             switch(method,
                    "WYht" = {
                      output <- WYht(X = X, f = f, theta = estSimpson, cmat = cmat, conf.level = conf.level, alternative = alternative, R = R, args = args)
                    },
                    "tsht" = {
                      output <- tsht(X = X, f = f, theta = estSimpsonf, cmat = cmat, conf.level = conf.level, alternative = alternative, R = R, args = args)
                    },
                    "rpht" = {
                      output <- rpht(X = X, f = f, theta = estSimpsonf, cmat = cmat, conf.level = conf.level, alternative = alternative, R = R, args = args)
                    },
                    "asht" = {
                      output <- asht(X = X, f = f, theta = estSimpsonf, cmat = cmat, conf.level = conf.level, alternative = alternative, args = args)
                    }
                    )
           }
           )
    return(output)
  }

