library(GeneralizedHyperbolic)
data(ghypParam)

test.ghypMean <- function(testParam = ghypSmallShape, n = 10000,
                          accuracy = 0.01) {
  for (i in 1:nrow(testParam)) {
    param <- testParam[i, ]

    # random number generation:
    x <- rghyp(n, param = param)

    # Compute mean of the sample data:
    sampleMean <- mean(x)

    # Get mean value from vgMean function:
    funMean <- ghypMean(param = param)

    # Precision within the accuracy value?
    difference <- abs(sampleMean - funMean)
    cat("param", sep = "\n")
    print(param)
    cat("sample mean", sep = "\n")
    print(sampleMean)
    cat("function mean", sep = "\n")
    print(funMean)
    cat("difference", sep = "\n")
    print(difference)
    cat("checkTrue", sep = "\n")
    print(checkTrue("test.ghypMean", difference < accuracy))
  }
}


test.ghypVar <- function(testParam = ghypSmallShape, n = 10000,
                         accuracy = 0.01) {
  for (i in 1:nrow(testParam)) {
    param <- testParam[i, ]

    # random number generation:
    x <- rghyp(n, param = param)

    # Compute variance of the sample data:
    sampleVar <- var(x)

    # Get mean value from vgVar function:
    funVar <- ghypVar(param = param)

    # Precision within the accuracy value?
    difference <- abs(sampleVar - funVar)
    cat("param", sep = "\n")
    print(param)
    cat("sample variance", sep = "\n")
    print(sampleVar)
    cat("function variance", sep = "\n")
    print(funVar)
    cat("difference", sep = "\n")
    print(difference)
    cat("checkTrue", sep = "\n")
    print(checkTrue("test.ghypVar", difference < accuracy))
  }
}

test.ghypSkew <- function(testParam = ghypSmallShape, n = 10000,
                          accuracy = 0.01) {
  for (i in 1:nrow(testParam)) {
    param <- testParam[i, ]

    # random number generation:
    x <- rghyp(n, param = param)

    # Compute skewness of the sample data:
    sampleSkew <- skewness(x)

    # Get skewness value from vgSkew function:
    funSkew <- ghypSkew(param = param)

    # Precision within the accuracy value?
    difference <- abs(sampleSkew - funSkew)
    cat("param", sep = "\n")
    print(param)
    cat("sample skewness", sep = "\n")
    print(sampleSkew)
    cat("function skewness", sep = "\n")
    print(funSkew)
    cat("difference", sep = "\n")
    print(difference)
    cat("checkTrue", sep = "\n")
    print(checkTrue("test.ghypSkew", difference < accuracy))
  }
}

test.ghypKurt <- function(testParam = ghypSmallShape, n = 10000,
                          accuracy = 0.01) {
  for (i in 1:nrow(testParam)) {
    param <- testParam[i, ]

    # random number generation:
    x <- rghyp(n, param = param)

    # Compute kurtosis of the sample data:
    sampleKurt <- kurtosis(x)

    # Get kurtosis value from vgKurt function:
    funKurt <- ghypKurt(param = param)

    # Precision within the accuracy value?
    difference <- abs(sampleKurt - funKurt)
    cat("param", sep = "\n")
    print(param)
    cat("sample kurtosis", sep = "\n")
    print(sampleKurt)
    cat("function kurtosis", sep = "\n")
    print(funKurt)
    cat("difference", sep = "\n")
    print(difference)
    cat("checkTrue", sep = "\n")
    print(checkTrue("test.ghypKurt", difference < accuracy))
  }
}

checkTrue <- function(f, bool) {
  if (bool) {
    paste(f, ": PASS", sep = "")
  } else {
    paste(f, ": FAIL", sep = "")
  }
}

test.ghypMean()
test.ghypVar()
test.ghypSkew()
test.ghypKurt()
