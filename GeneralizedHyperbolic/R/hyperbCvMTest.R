### Cramer-von Mises goodness of fit for the hyperbolic distribution
hyperbCvMTest <- function(x, mu = 0, delta = 1, alpha = 1, beta = 0,
                          param = c(mu, delta, alpha, beta),
                          conf.level = 0.95, ...) {

  if (!missing(conf.level) & (length(conf.level) != 1 |
                              !is.finite(conf.level) |
                              conf.level < 0 | conf.level > 1))
    stop("conf.level must be a single number between 0 and 1")

  if (length(param) != 4)
    stop("param vector must contain 4 values")

  DNAME <- deparse(substitute(x))
  NX <- length(x)

  if (NX < 2)
    stop("Not enough x observations")

  METHOD <- "Cramer-von Mises test of hyperbolic distribution"
  zvals <- phyperb(sort(x), param = param)
  STATISTIC <- sum((zvals - ((2 * (1:length(x)) - 1) / (2 * length(x))))^2) +
               1 / (12 * length(x))
  PARAMETER <- param
  names(STATISTIC) <- "Wsq"
  names(PARAMETER) <- c("mu", "delta", "alpha", "beta")
  names(METHOD) <- "Cramer-von Mises test of hyperbolic distribution"
  pValResult <- hyperbCvMTestPValue(delta, alpha, beta, STATISTIC)
  RVAL <- list(statistic = STATISTIC, method = METHOD, data.name = DNAME,
               parameter = PARAMETER, p.value = pValResult$pValue,
               warn = pValResult$warn)
  class(RVAL) <- "hyperbCvMTest"
  print(RVAL, ...)
} ## End hyperbCvMTest()

### Calculate P-Value of Cramer-von Mises test of the hyperbolic distribution
hyperbCvMTestPValue <- function(delta = 1, alpha = 1, beta = 0,
                                Wsq, digits = 3) {

  xi <- 1 / sqrt(1 + delta * sqrt(alpha^2 - beta^2))
  chi <- beta / (alpha * sqrt(1 + delta * sqrt(alpha^2 - beta^2)))

  xiList <- c(0.99, 0.95, seq(0.90, 0.1, by = -0.1))
  alphaList <- c(0.25, 0.1, 0.05, 0.025, 0.01)
  chiList <- seq(0, 0.8, by = 0.2)
  exactChi <- FALSE
  exactXi <- FALSE
  warn <- c(FALSE, FALSE)
  data(hyperbWSqTable)
  wsqTable <- hyperbWSqTable

  if (abs(chi) > xi)
    stop("Chi must be less than or equal to Xi")

  tol <- .Machine$double.eps

  if (chi > 0.8) {
    chi <- 0.79999
    warn[2] <- TRUE
  }

  for (j in 1:4) {
    if (chiList[j] < abs(chi) & (chiList [j + 1] - abs(chi)) > tol) {
      chiLo <- chiList[j]
      chiUp <- chiList[j + 1]
      wsqTable <- wsqTable[, j] + ((chi - chiList[j]) /
                                   (chiList[j + 1] - chiList[j])) *
                  (wsqTable[, j + 1] - wsqTable[, j])
    }
  }

  for (j in 1:5) {
    if (identical(all.equal(chiList[j], chi), TRUE)) {
      exactChi <- TRUE
      wsqTable <- wsqTable[, j]
    }
  }

  # correct xi values
  if (xi < 0.9 & abs(chi) > 0.6)
    xiList[4] <- xiList[4] + 0.01

  if (xi < 0.7 & abs(chi) > 0.4)
    xiList[6] <- xiList[6] + 0.01

  if (xi < 0.5 & abs(chi) > 0.2)
    xiList[8] <- xiList[8] + 0.01

  if (xi < 0.3 & abs(chi) > 0.0)
    xiList[10] <- xiList[10] + 0.01

  for (i in 1:11) {
    #print(i)
    if (xi > 0.99) {
      xi <- 0.99
      warn[1] <- TRUE
    }

    if (xi < 0.1) {
      xi <- 0.1
      warn[1] <- TRUE
    }

    if (identical(all.equal(xiList[i], xi), TRUE)) {
      exactXi <- TRUE
      wsqTable <- wsqTable[(5 * i - 4):(5 * i)]
    }

    if (xi < xiList[i] & xi - xiList[i + 1] > tol) {
      xiUp <- xiList[i]
      xiLo <- xiList[i+1]
      wsqTable <- wsqTable[(5 * i - 4):(5 * i + 5)]
    }
  }

  if (!exactXi) {
    wsqTable <- wsqTable[1:5] +
                ((xi - xiLo) / (xiUp - xiLo)) * (wsqTable[6:10] - wsqTable[1:5])
  }

  if (wsqTable[1] > Wsq)
    pValue <- "> 0.25"

  if (wsqTable[1] == Wsq)
    pValue <- "0.25"

  if (wsqTable[5] < Wsq)
    pValue <- "< 0.01"

  if (wsqTable[5] == Wsq)
    pValue <- "0.01"

  for (i in 1:4) {
    if (Wsq >= wsqTable[i] & Wsq < wsqTable[i + 1]) {
      pValue <- alphaList[i] + ((Wsq-wsqTable[i]) /
                                (wsqTable[i + 1] - wsqTable[i])) *
                               (alphaList[i + 1] - alphaList[i])
    }
  }

  if (is.numeric(pValue))
    pValue <- round(pValue, digits)

  pValResult <- list(pValue = pValue, warn = sum(warn))
  pValResult
} ## End hyperbCvMTestPValue()

### Print results of Cramer-von Mises goodness-of-fit test
### for the hyperbolic distribution
print.hyperbCvMTest <- function (x, prefix = "\t", ...) {

  cat("\n")
  writeLines(strwrap(x$method, prefix = "\t"))
  cat("\n")
  cat("data: ", x$data.name, "\n")
  out <- character()

  if (!is.null(x$statistic))
    out <- c(out, paste(names(x$statistic), "=",
             format(round(x$statistic, 4))))

  if (!is.null(x$parameter))
    out <- c(out, paste(names(x$parameter), "=",
             format(round(x$parameter, 3))))

  if (!is.null(x$p.value)) {
    fp <- as.character(x$p.value)
    out <- c(out, paste("p-value",
             if (substr(fp, 1, 1) == "<" | substr(fp, 1, 1) == ">") {
               fp
             } else {
               paste("=", fp)
             }))
  }

  if (x$warn != 0)
    warning("Estimated parameters are outside the table.\np-value may be incorrect")

  writeLines(strwrap(paste(out, collapse = ", ")))
  cat("\n")
  invisible(x)
} ## End of print.hyperbCvMTest()
