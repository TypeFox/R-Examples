#' Estimation of Nonlinear Regression Parameters with CRS4HC
#'
#'
#' @description
#' This function estimates the regression coefficients of a nonlinear regression function using least squares.
#' The minimization is performed by the CRS algorithm with four competing local heuristics. Algorithm is described in Tvrdík et al. (2007).
#'
#' @param formula (obligatory) a nonlinear \link{formula} including variables and parameters
#' @param data (obligatory) data frame in which to evaluate the variables in \code{formula}
#' @param a (obligatory) a vector of length equal to number of parameters representing lower bounds of search space (bounds for parameters must be specified in the same order they appear on right-hand side of \code{formula})
#' @param b (obligatory) a vector of length equal to number of parameters representing upper bounds of search space (bounds for parameters must be specified in the same order they appear on right-hand side of \code{formula})
#' @param N (optional) size of population. Default value is \code{10*length(a)}.
#' @param my_eps (optional) is used for stopping condition. Default value is 1e-15.
#' @param max_evals (optional) is used for stopping condition, specifies maximum number of objective function evaluations per dimension (dimension=nonlinear model parameter). Default value is 40000.
#' @param delta (optional) controls the competition of local heuristics. Default value is 0.05. delta > 0.
#' @param w0 (optional) controls the competition of local heuristics. Default value is 0.5. w0 > 0.
#' @usage crs4hc(formula, data, a, b, N, my_eps, max_evals, delta, w0)
#' @details There are implemented methods for generic functions \link{print}, \link{summary}, \link{plot}.
#'
#' @return
#' An S3 object of class \code{crs4hc}. This object is a list of:
#'   \item{model}{a list of two items, includes estimates of nonlinear model parameters and minimal residual sum of squares}
#'   \item{algorithmInfo}{a list of three items with some internal info about algorithm run}
#'   \item{data}{a data frame that was passed to function as the \code{data} argument}
#'   \item{other}{a list of four items which include info about nonlinear model \code{formula}}
#'
#' @export
#'
#' @importFrom stats runif as.formula
#'
#' @references
#'  Tvrdík, J., Křivý, I., and Mišík, L. Adaptive Population-based search:
#' Application to Estimation of Nonlinear Regression Parameters. \emph{Computational
#' Statistics and Data Analysis 52} (2007), 713–724. Preprint URL \url{http://www1.osu.cz/~tvrdik/wp-content/uploads/CSDA-06SAS03e.pdf}
#'
#' @examples
#' x <- c(1,2,3,5,7,10)
#' y <- c(109,149,149,191,213,224)
#' df <- data.frame(x=x, y=y)
#' lowerBounds <- c(1, 0.1)
#' upperBounds <- c(1000, 2)
#' mod <- crs4hc(y ~ b1 * (1-exp(-b2*x)), df, lowerBounds, upperBounds)
#' mod
#'
crs4hc <-
  function(formula,
           data ,
           a,
           b,
           N,
           my_eps = 1e-15,
           max_evals = 40000,
           delta = 0.05,
           w0 = 0.5)
  {
    if (missing(formula) || missing(data) || missing(a) || missing(b))
    {
      stop("One or more obligatory parameters are missing.")
    }

    dim <- length(a)
    noh <- 4

    if (dim != length(b))
    {
      stop("Length of vectors a and b are not equal.")
    }

    if (missing(N))
    {
      N <- 10 * dim
    }

    # formula handling
    formula <- as.formula(formula)
    if (!is.list(data) && !is.environment(data))
      stop("'data' must be a list or an environment")
    cll <- match.call()
    if (length(formula) == 2L) {
      formula[[3L]] <- formula[[2L]]
      formula[[2L]] <- 0
    }
    LHSVarNames <- all.vars(formula[[2]])
    if (LHSVarNames != formula[[2]])
    {
      data[, as.character(LHSVarNames)] <- eval(formula[[2]], envir = data)
      resVarName <- LHSVarNames
    }
    else
    {
      resVarName <- formula[[2]]
    }
    form2 <- formula
    form2[[2L]] <- 0
    RHS <- form2[[3]]
    varNamesRHS <- all.vars(form2)
    pnames <- varNamesRHS[is.na(match(varNamesRHS, colnames(data)))]

    remove <- c("pi", "e")
    pnames <- pnames[!pnames %in% remove]

    refl1count <- 0
    refl25count <- 0
    reflbcount <- 0
    cdeadpcount <- 0


    # Initialization of population matrix
    P <- matrix(data = 0,
                nrow = N,
                ncol = dim + 1)

    # 0th generation initialization
    for (i in 1:N)
    {
      P[i, (1:dim)] <- a + (b - a) * runif(dim)
      P[i, dim + 1] <- rss(RHS, P[i, 1:dim], data, resVarName, pnames)
    }

    ind_min <- which.min(P[, dim + 1])
    ind_max <- which.max(P[, dim + 1])
    fmin <- P[ind_min, dim + 1]
    fmax <- P[ind_max, dim + 1]

    num_of_evals <- N
    nrst <- 0
    wi <- vector(mode = "numeric", length = noh) + w0
    tss_part <-
      data[, as.character(resVarName)] - mean(data[, as.character(resVarName)])
    tss <- sum(tss_part ^ 2)
    tss_myeps <- tss * my_eps

    while ((fmax - fmin > tss_myeps) &&
           (num_of_evals < dim * max_evals))
    {
      hh = roulette_simple(wi)

      switch(hh,
        y <- refl_rwd(P, N, dim, c(2, 0.5)),
        y <- refl_rwd(P, N, dim, c(5, 1.5)) ,
        y <- refl_bestrwd(P, N, dim, 2, 0.5, ind_min, P[ind_min, 1:dim]),
        y <- cdeadp(P, N, dim, c(0.4, 0.9, fmin, fmax))
      )
      switch(hh,
        refl1count <- refl1count + 1,
        refl25count <- refl25count + 1,
        reflbcount <- reflbcount + 1,
        cdeadpcount <- cdeadpcount + 1
      )

      y <- zrcad(y, a, b)

      fy <- rss(RHS, y, data, resVarName, pnames)
      num_of_evals <- num_of_evals + 1


      if (fy < fmax)
      {
        P[ind_max, ] <- c(y, fy)
        if (fmax - fmin <= my_eps)
        {
          w <- w0
        }
        else
        {
          w <- (fmax - max(fy, fmin)) / (fmax - fmin)
        }
        if (is.nan(w))
        {
          w <- w0
        }
        wi[hh] <- wi[hh] + w
        p_min <- min(wi) / sum(wi)
        if (p_min < delta)
        {
          wi <- rep(0, noh) + w0
          nrst <- nrst + 1
        }
        ind_min <- which.min(P[, dim + 1])
        ind_max <- which.max(P[, dim + 1])
        fmin <- P[ind_min, dim + 1]
        fmax <- P[ind_max, dim + 1]
      }
    }
    b_star <- P[ind_min, 1:dim]
    names(b_star) <- pnames
    rss_star <- fmin

    result <- list()
    model <- list(estimates = b_star, rss = rss_star)
    result$model <- model
    result$call <- cll
    refl1Uses <- round((refl1count / num_of_evals) * 100, digits = 1)
    refl25Uses <- round((refl25count / num_of_evals) * 100, 1)
    reflbUses <- round((reflbcount / num_of_evals) * 100, 1)
    cdeadpUses <- round((cdeadpcount / num_of_evals) * 100, 1)
    result$algorithmInfo$heurUses <-
      c(
        REFL1 = refl1Uses,
        REFL25 = refl25Uses,
        REFLB = reflbUses,
        CDEADP = cdeadpUses
      )
    result$algorithmInfo$numOfEvals <- num_of_evals
    result$algorithmInfo$numOfResets <- nrst
    result$data <- data
    result$other$RHS <- RHS
    result$other$varNamesRHS <-
      varNamesRHS[!is.na(match(varNamesRHS, colnames(data)))]
    result$other$depVarName <- resVarName
    result$other$formula <- formula
    class(result) <- "crs4hc"
    return(result)
  }



#' Estimation of Nonlinear Regression Parameters with CRS4HCe
#'
#' @description
#' This function estimates the regression coefficients of a nonlinear regression function using least squares.
#' The minimization is performed by the CRS algorithm with four competing local heuristics and adaptive stopping condition. Algorithm is described in Tvrdík et al. (2007).
#'
#' @param formula (obligatory) a nonlinear \link{formula} including variables and parameters
#' @param data (obligatory) data frame in which to evaluate the variables in \code{formula}
#' @param a (obligatory) a vector of length equal to number of parameters representing lower bounds of search space (bounds for parameters must be specified in the same order they appear on right-hand side of \code{formula})
#' @param b (obligatory) a vector of length equal to number of parameters representing upper bounds of search space (bounds for parameters must be specified in the same order they appear on right-hand side of \code{formula})
#' @param N (optional) size of population. Default value is \code{10*length(a)}.
#' @param my_eps0 (optional) is used for adaptation of stopping condition. Default value is 1e-9.
#' @param gamma (optional) is used for adaptation of stopping condition. Default value is 1e7.
#' @param max_evals (optional) is used for stopping condition, specifies maximum number of objective function evaluations per dimension (dimension=nonlinear model parameter). Default values is 40000.
#' @param delta (optional) controls the competition of local heuristics. Default value is 0.05. delta > 0.
#' @param w0 (optional) controls the competition of local heuristics. Default value is 0.5. w0 > 0.
#' @usage crs4hce(formula, data , a, b, N, my_eps0, gamma, max_evals, delta, w0)
#'
#' @details It´s recommended to modify values of \code{my_eps0} and \code{gamma} together. There are implemented methods for generic functions \link{print}, \link{summary}, \link{plot}.
#'
#' @return
#' An S3 object of class \code{crs4hc}. This object is a list of:
#'   \item{model}{a list of two items, includes estimates of nonlinear model parameters and minimal residual sum of squares}
#'   \item{algorithmInfo}{a list of three items with some internal info about algorithm run}
#'   \item{data}{a data frame that was passed to function as the \code{data} argument}
#'   \item{other}{a list of four items which include info about nonlinear model \code{formula}}
#'
#'
#' @export
#'
#' @importFrom stats runif as.formula
#'
#' @references
#' Tvrdík, J., Křivý, I., and Mišík, L. Adaptive Population-based search:
#' Application to Estimation of Nonlinear Regression Parameters. \emph{Computational
#' Statistics and Data Analysis 52} (2007), 713–724. Preprint URL \url{http://www1.osu.cz/~tvrdik/wp-content/uploads/CSDA-06SAS03e.pdf}
#'
#' @examples
#' x <- c(1,2,3,5,7,10)
#' y <- c(109,149,149,191,213,224)
#' df <- data.frame(x=x, y=y)
#' lowerBounds <- c(1, 0.1)
#' upperBounds <- c(1000, 2)
#' mod <- crs4hce(y ~ b1 * (1-exp(-b2*x)), df, lowerBounds, upperBounds)
#' mod
#'
crs4hce <-
  function(formula,
           data ,
           a,
           b,
           N,
           my_eps0 = 1e-9,
           gamma = 1e7,
           max_evals = 40000,
           delta = 0.05,
           w0 = 0.5)
  {
    if (missing(formula) || missing(data) || missing(a) || missing(b))
    {
      stop("One or more obligatory parameters are missing.")
    }

    dim <- length(a)
    noh <- 4

    if (dim != length(b))
    {
      stop("Length of vectors a and b are not equal.")
    }

    if (missing(N))
    {
      N <- 10 * dim
    }


    formula <- as.formula(formula)
    if (!is.list(data) && !is.environment(data))
      stop("'data' must be a list or an environment")
    cll <- match.call()
    if (length(formula) == 2L) {
      formula[[3L]] <- formula[[2L]]
      formula[[2L]] <- 0
    }
    LHSVarNames <- all.vars(formula[[2]])
    if (LHSVarNames != formula[[2]])
    {
      data[, as.character(LHSVarNames)] <- eval(formula[[2]], envir = data)
      resVarName <- LHSVarNames
    }
    else
    {
      resVarName <- formula[[2]]
    }
    form2 <- formula
    form2[[2L]] <- 0
    RHS <- form2[[3]]
    varNamesRHS <- all.vars(form2)
    pnames <- varNamesRHS[is.na(match(varNamesRHS, colnames(data)))]

    remove <- c("pi", "e")
    pnames <- pnames[!pnames %in% remove]

    refl1count <- 0
    refl25count <- 0
    reflbcount <- 0
    cdeadpcount <- 0


    # Initialization of population matrix
    P <- matrix(data = 0,
                nrow = N,
                ncol = dim + 1)

    # 0th generation initialization
    for (i in 1:N)
    {
      P[i, (1:dim)] <- a + (b - a) * runif(dim)
      P[i, dim + 1] <- rss(RHS, P[i, 1:dim], data, resVarName, pnames)
    }

    ind_min <- which.min(P[, dim + 1])
    ind_max <- which.max(P[, dim + 1])
    fmin <- P[ind_min, dim + 1]
    fmax <- P[ind_max, dim + 1]

    num_of_evals <- N
    nrst <- 0
    wi <- vector(mode = "numeric", length = noh) + w0
    tss_part <-
      data[, as.character(resVarName)] - mean(data[, as.character(resVarName)])
    my_eps <- my_eps0
    tss <- sum(tss_part ^ 2)
    tss_myeps <- tss * my_eps
    one_minus_R2 <- 0
    pass <- F

    while (one_minus_R2 < (my_eps * gamma) && pass == F && gamma >= 10)
    {
      while ((fmax - fmin > tss_myeps) && (num_of_evals < dim * max_evals))
      {
        hh = roulette_simple(wi)

        switch(hh,
          y <- refl_rwd(P, N, dim, c(2, 0.5)),
          y <- refl_rwd(P, N, dim, c(5, 1.5)) ,
          y <-
            refl_bestrwd(P, N, dim, 2, 0.5, ind_min, P[ind_min, 1:dim]),
          y <- cdeadp(P, N, dim, c(0.4, 0.9, fmin, fmax))
        )
        switch(hh,
          refl1count <- refl1count + 1,
          refl25count <- refl25count + 1,
          reflbcount <- reflbcount + 1,
          cdeadpcount <- cdeadpcount + 1
        )

        y <- zrcad(y, a, b)

        fy <- rss(RHS, y, data, resVarName, pnames)
        num_of_evals <- num_of_evals + 1


        if (fy < fmax)
        {
          P[ind_max, ] <- c(y, fy)
          if (fmax - fmin <= my_eps)
          {
            w <- w0
          }
          else
          {
            w <- (fmax - max(fy, fmin)) / (fmax - fmin)
          }
          if (is.nan(w))
          {
            w <- w0
          }
          wi[hh] <- wi[hh] + w
          p_min <- min(wi) / sum(wi)
          if (p_min < delta)
          {
            wi <- rep(0, noh) + w0
            nrst <- nrst + 1
          }
          ind_min <- which.min(P[, dim + 1])
          ind_max <- which.max(P[, dim + 1])
          fmin <- P[ind_min, dim + 1]
          fmax <- P[ind_max, dim + 1]
        }
        pass <- T
      }
      if (pass == F)
      {
        gamma <- gamma / 10
      }
      one_minus_R2 <- fmin / tss

      if (one_minus_R2 < (my_eps0 * gamma) && pass == T)
      {
        my_eps <- my_eps * 0.1
        tss_myeps <- tss * my_eps
        pass <- F
      }
    }
    b_star <- P[ind_min, 1:dim]
    if (pass)
    {
      my_eps <- my_eps * 10
    }
    names(b_star) <- pnames
    rss_star <- fmin

    result <- list()
    model <- list(estimates = b_star, rss = rss_star)
    result$model <- model
    result$call <- cll
    refl1Uses <- round((refl1count / num_of_evals) * 100, digits = 1)
    refl25Uses <- round((refl25count / num_of_evals) * 100, 1)
    reflbUses <- round((reflbcount / num_of_evals) * 100, 1)
    cdeadpUses <- round((cdeadpcount / num_of_evals) * 100, 1)
    result$algorithmInfo$heurUses <-
      c(
        REFL1 = refl1Uses,
        REFL25 = refl25Uses,
        REFLB = reflbUses,
        CDEADP = cdeadpUses
      )
    result$algorithmInfo$numOfEvals <- num_of_evals
    result$algorithmInfo$numOfResets <- nrst
    result$algorithmInfo$my_eps <- my_eps
    result$data <- data
    result$other$RHS <- RHS
    result$other$varNamesRHS <-
      varNamesRHS[!is.na(match(varNamesRHS, colnames(data)))]
    result$other$depVarName <- resVarName
    result$other$formula <- formula
    class(result) <- "crs4hc"
    return(result)
  }

#' @importFrom stats runif
cdeadp <- function(P, N, d, vecpar)
{
  Fmin <- vecpar[1]
  CR <- vecpar[2]
  fmin <- vecpar[3]
  fmax <- vecpar[4]
  vyb <- nahvyb(N, 4)
  r1 <- P[vyb[1], 1:d]
  r2 <- P[vyb[2], 1:d]
  r3 <- P[vyb[3], 1:d]
  y <- P[vyb[4], 1:d]
  pom1 <- Fmin
  pom2 <- 1
  pom3 <- 1

  if (abs(fmin) > 0)
  {
    pom2 <- abs(fmax / fmin)
  }
  if (pom2 < 1)
  {
    pom1 <- 1 - pom2
  }
  else if (abs(fmax) > 0)
  {
    pom1 <- 1 - abs(fmin / fmax)
  }
  ff <- max(Fmin, pom1)
  v <- r1 + ff * (r2 - r3)
  change <- which(runif(d) < CR)
  if (length(change) == 0)
  {
    change <- 1 + trunc(d * runif(1))
  }
  y[change] <- v[change]
  return(y)
}

#' @importFrom stats runif
nahvyb <- function(N, k)
{
  opora <- 1:N
  nahv <- rep(0, k)

  for (i in 1:k)
  {
    index <- 1 + trunc(runif(1) * length(opora))
    nahv[i] <- opora[index]
    opora <- opora[-index]
  }
  return(nahv)
}

#' @importFrom stats runif
nahvyb_expt <- function(N, k, expt)
{
  opora <- 1:N
  if (nargs() == 3)
  {
    opora <- opora[-expt]
  }
  vyb <- rep(0, k)

  for (i in 1:k)
  {
    index <- 1 + trunc(runif(1) * length(opora))
    vyb[i] <- opora[index]
    opora <- opora[-index]
  }
  return(vyb)
}

#' @importFrom stats runif
refl_bestrwd <- function(P, N, d, alpha, shift_d, indmin, xmin)
{
  vyb <- nahvyb_expt(N, d, indmin)
  S <- P[vyb, ]
  indx <- which.max(S[, d + 1])
  x <- S[indx, 1:d]
  S <- S[-indx, ]
  cs <- vector(mode = "numeric", length = 0)
  for (i in 1:d)
  {
    cs <- c(cs, sum(S[i]))
  }
  g <- (xmin + cs) / d
  y <- g + (g - x) * (shift_d + (alpha - 2 * shift_d) * runif(1))
  return(y)
}

#' @importFrom stats runif
refl_rwd <- function(P, N, d, alpha)
{
  shift_d <- alpha[2]
  alpha <- alpha[1]
  vyb <- nahvyb(N, d + 1)
  S <- P[vyb, ]
  indx <- which.max(S[, d + 1])
  x <- S[indx, 1:d]
  S <- S[-indx, ]
  g <- colMeans(S[, 1:d])
  y <- g + (g - x) * (shift_d + (alpha - 2 * shift_d) * runif(1))
  return(y)
}

#' @importFrom stats runif
roulette_simple <- function(cutpoints)
{
  h <- length(cutpoints)
  ss <- sum(cutpoints)
  cp <- vector(mode = "numeric", length = h)
  cp[1] <- cutpoints[1]

  for (i in 2:h)
  {
    cp[i] <- cp[i - 1] + cutpoints[i]
  }
  cp <- cp / ss
  res <- 1 + trunc(sum(cp < runif(1)))
}

rss <- function(formula, betas, data, resVarName, pnames)
{
  resVarName <- as.character(resVarName)
  names(betas) <- pnames
  ls <- list()
  ls <- append(ls, as.list(data[, 1:length(data[1, ])]))
  ls <- append(ls, as.list(betas))

  y_hat <- eval(formula, envir = ls)
  result <- sum((data[, resVarName] - y_hat) ^ 2)
  return(result)
}

zrcad <- function(y, a, b)
{
  zrc <- c(which(y < a), which(y > b))
  for (i in zrc)
  {
    while ((y[i] < a[i]) || (y[i] > b[i]))
    {
      if (y[i] > b[i])
        y[i] <- 2 * b[i] - y[i]
      else
      {
        if (y[i] < a[i])
          y[i] <- 2 * a[i] - y[i]
      }
    }
  }
  return(y)
}



#' @export
print.crs4hc <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)
  cat("Parameters estimates:\n")
  print(x$mode$estimates)
  cat("Residual Sum of Squares:\n")
  print(x$model$rss)
}


#' @export
summary.crs4hc  <- function(object, ...)
{
  cat(
    "********************************************************************************\n"
  )
  cat("MODEL\n")
  print.crs4hc(object)
  cat(
    "********************************************************************************\n"
  )
  cat("INTERNAL INFO\n")
  cat("Heuristics relative uses: \n")
  print(object$algorithmInfo$heurUses)
  cat("Number of objective function evaluation:\n")
  print(object$algorithmInfo$numOfEvals)
  cat("Number of probability resets:\n")
  print(object$algorithmInfo$numOfResets)
}


#' @export
#' @importFrom graphics par plot points
plot.crs4hc <- function(x, ...)
{
  depVarName <- x$other$depVarName
  numOfIndVars <- length(x$other$varNamesRHS)
  if (numOfIndVars > 1)
  {
    par(mfcol = c(numOfIndVars / 2, 2))
  }
  else
  {
    par(mfrow = c(1, 1))
  }

  for (i in 1:numOfIndVars)
  {
    indVarName <- x$other$varNamesRHS[i]
    x_min <- min(x$data[indVarName])
    x_max <- max(x$data[indVarName])
    x_seq <- seq(from = x_min,
                 to = x_max,
                 by = (x_max - x_min) / 100)
    if (i == 1)
    {
      df <- data.frame(x = x_seq)
    }
    else
    {
      df <- cbind(df, x = x_seq)
    }
    colnames(df) <- x$other$varNamesRHS[1:i]
  }
  ls <- list()
  ls <- append(ls, as.list(df))
  ls <- append(ls, as.list(x$model$estimates))
  evals <- eval(x$other$RHS, ls)
  df[, "y"] <- evals
  colnames(df) <- c(x$other$varNamesRHS, depVarName)
  for (i in 1:numOfIndVars)
  {
    indVarName <- x$other$varNamesRHS[i]
    plot(df[, as.character(indVarName)],
         df[, as.character(depVarName)],
         type = "l",
         xlab = indVarName,
         ylab = depVarName)
    points(x = x$data[, as.character(indVarName)],
           y = x$data[, as.character(depVarName)],
           col = "red")
  }
  par(mfrow = c(1, 1))
}
